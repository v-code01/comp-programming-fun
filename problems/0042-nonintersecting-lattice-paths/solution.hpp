#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0042 {

using i64 = std::int64_t;

constexpr i64 kMod = 1000000007;  // prime -- every nonzero residue has an inverse.

// coordinates run to 2e5, so the longest monotone path is (2e5)+(2e5) steps and
// the biggest binomial we ever ask for is C(4e5, ·). size the tables once, for
// the worst case, and never think about it again.
constexpr int kMaxCoord = 200000;
constexpr int kMaxTop = 2 * kMaxCoord;  // 4e5.

struct Pt {
    int x, y;
};

inline bool operator==(Pt a, Pt b) { return a.x == b.x && a.y == b.y; }

// (a*b) mod p. INVARIANT, and every caller here holds it: a and b are already
// reduced, so both are in [0, p) with p < 2^30. the product is then under 2^60 and
// lands in an int64 exactly -- the static_assert is the proof, not a hope.
//
// the __int128 version of this is the reflex, and it is 4x slower: a 128-bit
// modulo can't use the compiler's constant-divisor magic-number trick, so it
// becomes a __umodti3 libcall. 4.26 ns against 1.08 ns, measured. this is the hot
// instruction of the whole solve -- k^3/3 of them in the elimination -- so the
// int64 form is not a micro-optimization, it's the difference.
static_assert((kMod - 1) <= 3037000499LL, "p^2 must fit in int64 for mul() to be exact");

inline i64 mul(i64 a, i64 b) { return a * b % kMod; }

// square-and-multiply. only ever called with e = p-2, so ~30 iterations.
inline i64 pw(i64 a, i64 e) {
    a %= kMod;
    if (a < 0) a += kMod;
    i64 r = 1;
    while (e) {
        if (e & 1) r = mul(r, a);
        a = mul(a, a);
        e >>= 1;
    }
    return r;
}

// fermat: p is prime, so a^(p-2) == a^{-1} (mod p). this is the line that lets
// gaussian elimination divide by a pivot. over a composite modulus the whole
// determinant routine below falls apart -- the pivot might not be invertible and
// the field argument evaporates.
inline i64 inv(i64 a) { return pw(a, kMod - 2); }

// factorials and their inverses to 4e5. built once, on first use, O(maxCoord).
// the inverse table walks down from a single fermat inverse -- fi[i-1] = fi[i]*i
// -- so we pay one modpow, not 4e5 of them.
struct Tables {
    std::vector<i64> f, fi;

    Tables() : f(kMaxTop + 1), fi(kMaxTop + 1) {
        f[0] = 1;
        for (int i = 1; i <= kMaxTop; ++i) f[i] = mul(f[i - 1], i);
        fi[kMaxTop] = inv(f[kMaxTop]);
        for (int i = kMaxTop; i >= 1; --i) fi[i - 1] = mul(fi[i], i);
    }
};

inline const Tables& tables() {
    static const Tables t;
    return t;
}

// pay the O(maxCoord) table build up front, so a caller timing a solve is timing
// the solve.
inline void prepare() { (void)tables(); }

inline i64 binom(int n, int r) {
    if (r < 0 || n < 0 || r > n) return 0;
    const Tables& t = tables();
    return mul(t.f[n], mul(t.fi[r], t.fi[n - r]));
}

// monotone paths from a to b, moving only right and up. you take dx rights and
// dy ups in some order, so the count is C(dx+dy, dx) -- choose which of the
// steps are the rights.
//
// b strictly left of a, or strictly below it, and there is no path at all. that
// zero is load-bearing: it is what makes the matrix triangular-ish in the tight
// configurations and what drives the answer to 0 when an end is unreachable.
inline i64 path_count(Pt a, Pt b) {
    if (b.x < a.x || b.y < a.y) return 0;
    const int dx = b.x - a.x;
    const int dy = b.y - a.y;
    return binom(dx + dy, dx);
}

// determinant over GF(p) by gaussian elimination. O(k^3).
//
// row reduce to upper triangular, multiply the pivots. two things to get right,
// and they are the two things that break in every buggy version of this:
//   1. the pivot has to be inverted, not divided by. fermat gives the inverse.
//   2. every row swap flips the sign of the determinant. track it -- forget it
//      and you get the right magnitude with the wrong sign, which mod p looks
//      like a completely unrelated number.
// no pivot in a column means the column is dependent on the ones before it.
// determinant 0, stop early.
inline i64 det_mod(std::vector<std::vector<i64>> m) {
    const int n = static_cast<int>(m.size());
    i64 res = 1;

    for (int c = 0; c < n; ++c) {
        int piv = -1;
        for (int r = c; r < n; ++r) {
            if (m[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] != 0) {
                piv = r;
                break;
            }
        }
        if (piv < 0) return 0;  // singular -- the column is a combination of the left ones.

        if (piv != c) {
            std::swap(m[static_cast<std::size_t>(piv)], m[static_cast<std::size_t>(c)]);
            res = (kMod - res) % kMod;  // one transposition, one sign flip.
        }

        // hoist the pivot row. m is a vector of vectors, so leaving it indexed in
        // the inner loop means chasing a pointer per cell for a row that never
        // moves.
        const std::vector<i64>& prow = m[static_cast<std::size_t>(c)];
        const i64 p = prow[static_cast<std::size_t>(c)];
        res = mul(res, p);
        const i64 ip = inv(p);  // one modpow per column, k of them total. noise next to k^3.

        for (int r = c + 1; r < n; ++r) {
            std::vector<i64>& row = m[static_cast<std::size_t>(r)];
            const i64 factor = row[static_cast<std::size_t>(c)];
            if (factor == 0) continue;  // already eliminated. the tight configs are full of these.
            const i64 k = mul(factor, ip);
            for (int j = c; j < n; ++j) {
                // mul() already reduced its product, so cell - t sits in
                // (-p, p). one conditional add lands it back in [0, p) -- no second
                // 64-bit division. that division was half the inner loop.
                const i64 t = row[static_cast<std::size_t>(j)] -
                              mul(k, prow[static_cast<std::size_t>(j)]);
                row[static_cast<std::size_t>(j)] = t < 0 ? t + kMod : t;
            }
        }
    }
    return res;
}

// --- lindstrom-gessel-viennot ---
//
// build M[i][j] = number of monotone paths a_i -> b_j. the answer is det(M).
//
// WHY. expand the determinant by its definition:
//
//     det(M) = sum over permutations sigma of  sign(sigma) * prod_i M[i][sigma(i)]
//            = sum over sigma  sign(sigma) * (families of k paths, path i running
//                                             from a_i to b_sigma(i))
//
// that inner product is literally a count of path FAMILIES -- every family, the
// tangled ones included. so the determinant is a signed sum over all families,
// under all pairings. the lemma says all the tangled ones kill each other.
//
// the involution. take any family that has an intersection. fix the smallest i
// whose path P_i touches some other path. walk P_i from a_i and stop at the first
// lattice point it shares with another path; among the paths through that point
// take the smallest index j > i. now swap the TAILS of P_i and P_j at that point.
//
//   - the weight is unchanged. the same multiset of edges is still on the board;
//     you only relabeled which path owns which tail.
//   - the pairing changed by exactly one transposition: i now ends at
//     b_sigma(j) and j at b_sigma(i). so sign(sigma) flips.
//   - applying it twice gets you back. the choice of i, of the point, and of j is
//     made from data the swap does not disturb (the prefixes of P_i and P_j up to
//     that point are untouched, and so is every other path), so it is an
//     involution on the intersecting families.
//
// an involution that flips sign and preserves weight cancels everything it
// touches, in pairs. what survives the sum is exactly the families with NO
// intersection.
//
// COMPATIBILITY -- the condition that turns the signed sum into a plain count.
// the surviving non-intersecting families still carry sign(sigma), so they'd only
// add up to a count if every one of them has sigma = identity. that is the
// standing assumption on the input, and it is not free:
//
//     a_1..a_k sorted with x nondecreasing and y nonincreasing,
//     b_1..b_k sorted the same way.
//
// each list walks south-east. under that ordering, any non-identity pairing has
// an inversion i < l with sigma(i) > sigma(l), which routes a path from the more
// north-west start to the more south-east end while another goes from the more
// south-east start to the more north-west end. two monotone staircases in that
// arrangement cannot miss each other -- and on the integer lattice, two monotone
// staircases that cross have to meet AT a lattice point (axis-aligned unit
// segments on integer lines only ever intersect at integer points). so every
// sigma != id contributes zero non-intersecting families, and
//
//     det(M) = number of vertex-disjoint families with path i from a_i to b_i.
//
// that is the whole lemma. drop the ordering and det(M) is still a perfectly good
// signed count -- it just stops being THE count.
//
// LOWER BOUND, honestly. you have to read the 2k points, so nothing beats
// Omega(k). we spend O(k^2) binomials to fill M and O(k^3) to reduce it, plus a
// one-time O(maxCoord) for the factorial tables. the alternative that answers the
// question directly -- enumerate path families and throw out the tangled ones --
// is exponential in the path lengths (that's the oracle in reference.hpp, and it
// dies around span 4). LGV is the trade: an exponential count collapses into a
// cubic one, because the cancellation does the enumeration for you.
inline i64 solve(const std::vector<Pt>& a, const std::vector<Pt>& b) {
    const int k = static_cast<int>(a.size());
    std::vector<std::vector<i64>> m(static_cast<std::size_t>(k),
                                    std::vector<i64>(static_cast<std::size_t>(k), 0));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                path_count(a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(j)]);
    return det_mod(std::move(m));
}

}  // namespace p0042
