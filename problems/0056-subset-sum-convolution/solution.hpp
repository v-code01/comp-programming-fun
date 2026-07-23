#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0056 {

// the field. 998244353 is prime, but we never divide here -- only add, subtract,
// and multiply. the modulus is just the ring the coefficients live in.
inline constexpr std::uint32_t kMod = 998244353;

// --- subset sum convolution ---
//
// c[S] = sum over T subset of S of a[T] * b[S\T]   (mod p).
// this is a DISJOINT-union convolution: T and S\T never share a bit, and their
// union is exactly S. the naive read is "for each S, walk its submasks" -- that
// is O(3^k) (each bit is in T, in S\T, or outside S), and at k=20 that is ~3.5e9
// with a multiply inside. dead.
//
// the trick is to make the disjointness visible to a plain subset-sum transform,
// which on its own cannot see it. split everything by POPCOUNT (rank).
//
// step 1 -- rank-layered zeta.
//   fa_j[T] = a[T] if popcount(T)==j else 0. run the SOS/zeta transform on each
//   of the k+1 layers separately:
//       ahat_j[S] = sum over T subset of S of fa_j[T]
//                 = sum over T subset of S with popcount(T)==j of a[T].
//   same for bhat. this is (k+1) transforms, each O(2^k * k).
//
// step 2 -- rank convolution, pointwise per S.
//   for each S and each target rank i, take the degree-i slice:
//       chat_i[S] = sum over j=0..i of ahat_j[S] * bhat_{i-j}[S].
//   pointwise product of two zeta transforms is the zeta of the covering
//   product, so chat_i = zeta(D_i) where
//       D_i[V] = sum over T,U with T union U == V and popcount(T)+popcount(U)==i
//                of a[T] * b[U].
//   nothing forces T and U disjoint yet -- overlap is still allowed.
//
// step 3 -- mobius-invert each layer, then read the rank that matches S.
//   invert the SOS on layer chat_i to recover D_i. now the disjointness falls
//   out: for a mask V with popcount(V)==i,
//       |T union U| = |V| = i  and  |T|+|U| = i  =>  |T intersect U| = 0,
//   so T and U are disjoint and U = V\T. hence
//       D_i[V] = sum over T subset of V of a[T] * b[V\T] = c[V]   exactly.
//   for popcount(V) != i the overlap terms poison D_i[V], so we never use them.
//   c[S] = D_{popcount(S)}[S] -- pick the layer whose rank equals popcount(S).
//   that is the whole reason the rank index acts like the degree variable a plain
//   SOS is blind to: it recovers the disjointness the covering product loses.
//
// complexity. reading a,b is Omega(2^k); the output is 2^k values so Omega(2^k).
// the k+1 rank layers cost (k+1) zeta transforms at O(2^k * k) each, the rank
// convolution is O(2^k * k^2), and the k+1 mobius inversions are again O(2^k * k)
// each. total O(2^k * k^2), against O(3^k) for the naive submask walk. the extra
// k^2 is the price of restoring disjointness through the rank index.
//
// space O(2^k * k): the k+1 zeta layers of a, of b, and of the result. values are
// stored as uint32 (every residue is < 2^30), which halves the footprint and the
// bandwidth versus 64-bit.

namespace detail {

// subset-sum (zeta) transform in place over one rank layer:
//   f[S] <- sum over T subset of S of f[T].
// for each bit, fold the "bit clear" sibling into the "bit set" mask. each add is
// of two residues < p < 2^30, so the sum is < 2^31 and never overflows uint32 --
// one conditional subtract normalizes it, cheaper than a modulo.
inline void zeta(std::uint32_t* f, int k, std::size_t size) {
    for (int i = 0; i < k; ++i) {
        const std::size_t bit = static_cast<std::size_t>(1) << i;
        for (std::size_t mask = 0; mask < size; ++mask) {
            if (mask & bit) {
                const std::uint32_t v = f[mask] + f[mask ^ bit];
                f[mask] = v >= kMod ? v - kMod : v;
            }
        }
    }
}

// inverse subset-sum (mobius) transform in place, the exact undo of zeta:
//   f[S] <- sum over T subset of S of (-1)^(|S|-|T|) f[T].
// subtract the sibling instead of adding it. f[mask^bit] < p, so f[mask]+p never
// overflows uint32 and the difference stays non-negative before the normalize.
inline void mobius(std::uint32_t* f, int k, std::size_t size) {
    for (int i = 0; i < k; ++i) {
        const std::size_t bit = static_cast<std::size_t>(1) << i;
        for (std::size_t mask = 0; mask < size; ++mask) {
            if (mask & bit) {
                const std::uint32_t v = f[mask] + kMod - f[mask ^ bit];
                f[mask] = v >= kMod ? v - kMod : v;
            }
        }
    }
}

}  // namespace detail

// c[S] = sum over T subset of S of a[T]*b[S\T], mod p, for every S in [0, 2^k).
// a and b are indexed by subset, each entry already a residue in [0, p). k <= 20.
inline std::vector<std::uint32_t> solve(const std::vector<std::uint32_t>& a,
                                        const std::vector<std::uint32_t>& b,
                                        int k) {
    const std::size_t size = static_cast<std::size_t>(1) << k;
    const int layers = k + 1;  // ranks 0..k inclusive.

    // rank-major layout [rank*size + S]: each layer is contiguous, so the zeta and
    // mobius sweeps stay cache-friendly. the pointwise convolution pays a strided
    // gather per S, but it reads only 2*(k+1) cells there and works from locals.
    std::vector<std::uint32_t> ahat(static_cast<std::size_t>(layers) * size, 0);
    std::vector<std::uint32_t> bhat(static_cast<std::size_t>(layers) * size, 0);

    // scatter a,b into their popcount layer -- fa_j[S] = a[S] iff popcount(S)==j.
    for (std::size_t S = 0; S < size; ++S) {
        const int r = __builtin_popcountll(static_cast<unsigned long long>(S));
        ahat[static_cast<std::size_t>(r) * size + S] = a[S];
        bhat[static_cast<std::size_t>(r) * size + S] = b[S];
    }

    // step 1: zeta each rank layer of a and of b.
    for (int r = 0; r < layers; ++r) {
        detail::zeta(ahat.data() + static_cast<std::size_t>(r) * size, k, size);
        detail::zeta(bhat.data() + static_cast<std::size_t>(r) * size, k, size);
    }

    // step 2: rank convolution, pointwise per S.
    //   chat_i[S] = sum_{j=0..i} ahat_j[S] * bhat_{i-j}[S].
    // each product is < p^2 < 2^60 and we reduce it to < p before adding, so the
    // running sum of at most k+1 residues stays under (k+1)*p < 2^35 -- no uint64
    // overflow, one final reduce per (S, i).
    std::vector<std::uint32_t> chat(static_cast<std::size_t>(layers) * size, 0);
    std::vector<std::uint32_t> av(static_cast<std::size_t>(layers));
    std::vector<std::uint32_t> bv(static_cast<std::size_t>(layers));
    for (std::size_t S = 0; S < size; ++S) {
        for (int r = 0; r < layers; ++r) {
            av[static_cast<std::size_t>(r)] =
                ahat[static_cast<std::size_t>(r) * size + S];
            bv[static_cast<std::size_t>(r)] =
                bhat[static_cast<std::size_t>(r) * size + S];
        }
        for (int i = 0; i < layers; ++i) {
            std::uint64_t acc = 0;
            for (int j = 0; j <= i; ++j) {
                acc += static_cast<std::uint64_t>(av[static_cast<std::size_t>(j)]) *
                       bv[static_cast<std::size_t>(i - j)] % kMod;
            }
            chat[static_cast<std::size_t>(i) * size + S] =
                static_cast<std::uint32_t>(acc % kMod);
        }
    }

    // ahat/bhat are dead now -- drop them before the mobius pass so peak memory is
    // two layer-arrays, not three.
    ahat.clear();
    ahat.shrink_to_fit();
    bhat.clear();
    bhat.shrink_to_fit();

    // step 3: mobius-invert each layer, then read the layer whose rank == |S|.
    for (int r = 0; r < layers; ++r) {
        detail::mobius(chat.data() + static_cast<std::size_t>(r) * size, k, size);
    }

    std::vector<std::uint32_t> c(size);
    for (std::size_t S = 0; S < size; ++S) {
        const int r = __builtin_popcountll(static_cast<unsigned long long>(S));
        c[S] = chat[static_cast<std::size_t>(r) * size + S];
    }
    return c;
}

}  // namespace p0056
