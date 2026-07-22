#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0052 {

using i64 = std::int64_t;
using u64 = std::uint64_t;

// GF(p), p = 998244353 = 119 * 2^23 + 1. the 2^23 is the point -- p-1 holds a
// factor of 2^23, so a primitive 2^k-th root of unity lives in the field for
// every k up to 23. g = 3 generates the multiplicative group, so g^((p-1)/2^k)
// has exact order 2^k. that's the exact-arithmetic stand-in for e^{-2pi i / N}
// the complex FFT uses -- same butterflies, integer end to end, no rounding to
// argue about. same p and g as 0043, because the multiply, the divide, and the
// product tree all lean on the same NTT.
inline constexpr i64 kMod = 998244353;
inline constexpr i64 kG = 3;  // primitive root of GF(998244353).

// (a*b) mod p. both factors sit under p < 2^30, so the product tops near 2^60 --
// fits a u64 with room, no __int128. WATCH: callers must keep operands in [0, p);
// hand this a raw point and the u64 cast still holds but the residue lies. % by
// the constexpr kMod lowers to a magic-multiply, not a hardware divide -- the
// difference in a butterfly that runs tens of millions of times.
inline i64 mul(i64 a, i64 b) {
    return static_cast<i64>(static_cast<u64>(a) * static_cast<u64>(b) %
                            static_cast<u64>(kMod));
}

// a^e mod p, square-and-multiply. roots of unity and the fermat inverse ride it.
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

// fermat: p prime, so a^(p-2) == a^{-1} (mod p). the field is the whole reason
// division exists here -- newton divides by the leading coeff to seed, the NTT
// divides by the size to invert. over a composite modulus none of it holds.
inline i64 inv(i64 a) { return pw(a, kMod - 2); }

// iterative in-place radix-2 NTT. n a power of two, inputs already in [0, p).
// invert=false forward, invert=true inverse (same butterflies, inverse root,
// then a divide by n). three parts, same as the complex FFT: bit-reversal
// permutation, butterflies over powers of the len-th root, and the 1/n scale on
// the way back.
inline void ntt(std::vector<i64>& a, bool invert) {
    const int n = static_cast<int>(a.size());

    // bit-reversal permutation. j mirrors i by carrying the top-bit flip down --
    // no per-element ctz, just an incrementing reversed counter.
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // stage over block sizes len = 2, 4, …, n. each block folds its two halves.
    for (int len = 2; len <= n; len <<= 1) {
        // primitive len-th root. len divides 2^23 divides p-1, so (p-1)/len is
        // an integer and g^that has exact order len.
        i64 w = pw(kG, (kMod - 1) / len);
        if (invert) w = inv(w);  // inverse transform walks the roots backward.
        const int half = len >> 1;
        for (int i = 0; i < n; i += len) {
            i64 wn = 1;  // w^0, w^1, … across the half-block. exact, no drift.
            for (int j = 0; j < half; ++j) {
                const i64 u = a[i + j];
                const i64 v = mul(a[i + j + half], wn);
                i64 s = u + v;
                if (s >= kMod) s -= kMod;  // u,v < p -> u+v < 2p, one subtract.
                i64 d = u - v;
                if (d < 0) d += kMod;  // u-v in (-p, p), one add.
                a[i + j] = s;
                a[i + j + half] = d;
                wn = mul(wn, w);
            }
        }
    }

    if (invert) {
        const i64 ninv = inv(n);  // undo the factor of n the round introduced.
        for (auto& x : a) x = mul(x, ninv);
    }
}

// linear convolution over GF(p): res[k] = sum_i a[i] b[k-i] mod p. pad to a
// power of two big enough to hold the full product, transform, multiply
// pointwise, invert. exact -- integer arithmetic the whole way. this is the one
// primitive the product tree and the divide both spend their time in.
inline std::vector<i64> convolve(std::vector<i64> a, std::vector<i64> b) {
    if (a.empty() || b.empty()) return {};
    const int need = static_cast<int>(a.size() + b.size()) - 1;
    int sz = 1;
    while (sz < need) sz <<= 1;
    a.resize(static_cast<std::size_t>(sz), 0);
    b.resize(static_cast<std::size_t>(sz), 0);
    ntt(a, false);
    ntt(b, false);
    for (int i = 0; i < sz; ++i) a[i] = mul(a[i], b[i]);
    ntt(a, true);
    a.resize(static_cast<std::size_t>(need));
    return a;
}

// B(x) = A(x)^{-1} mod x^n: the unique degree-<n poly with A*B == 1 mod x^n.
// requires a_0 != 0. NEWTON DOUBLING -- hold B correct mod x^k, push to mod
// x^{2k} in one step: B' = B*(2 - A*B). WHY it doubles: with E = 1 - A*B == 0
// mod x^k, A*B' = (1-E)(1+E) = 1 - E^2, and E^2 == 0 mod x^{2k}. precision k ->
// 2k every step, seeded at k=1 by B = [a_0^{-1}]. O(n log n) total -- the k's
// double so the transform costs form a geometric series in the last term.
//
// here it is the engine of polynomial division: the reversed divisor is invertible
// (its constant term is the divisor's leading coeff), and one inverse turns long
// division into a convolution. see poly_mod.
inline std::vector<i64> inverse(const std::vector<i64>& A, int n) {
    i64 a0 = A[0] % kMod;
    if (a0 < 0) a0 += kMod;
    std::vector<i64> b{inv(a0)};

    for (int len = 1; len < n; len <<= 1) {
        const int t = len << 1;  // target precision this step: mod x^t.

        // A truncated to x^t -- terms at degree >= t never touch B mod x^t.
        std::vector<i64> at(A.begin(),
                            A.begin() + (static_cast<int>(A.size()) < t
                                             ? static_cast<int>(A.size())
                                             : t));

        // C = A*B, then D = 2 - C, both mod x^t.
        std::vector<i64> c = convolve(at, b);
        c.resize(static_cast<std::size_t>(t), 0);
        for (auto& x : c) x = x ? kMod - x : 0;  // D = -C …
        c[0] += 2;                                // … + 2.
        if (c[0] >= kMod) c[0] -= kMod;

        // B' = B*D mod x^t.
        std::vector<i64> nb = convolve(b, c);
        nb.resize(static_cast<std::size_t>(t), 0);
        b = std::move(nb);
    }

    b.resize(static_cast<std::size_t>(n), 0);
    return b;
}

// strip trailing (high-degree) zero coefficients so degree reads true.
inline void trim(std::vector<i64>& a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

// R = A mod B, deg R < deg B. B monic (b_top == 1), deg B >= 1 -- exactly the
// shape every product-tree node has. this is polynomial long division done fast.
//
// A = Q*B + R, deg Q = deg A - deg B, deg R < deg B. REVERSED-POLY TRICK: write
// A^R(x) = x^{deg A} A(1/x), the coefficient vector flipped. from A = QB + R,
// with deg R < deg B every term of x^{deg A} R(1/x) sits at degree > deg A - deg B,
// so modulo x^{deg A - deg B + 1} they vanish and
//
//   A^R == Q^R * B^R   (mod x^{dq+1}),  dq = deg A - deg B.
//
// B^R has constant term b_top -- monic here, so it's 1, hence invertible. one
// newton inverse gives Q^R = A^R * (B^R)^{-1} mod x^{dq+1}; flip to get Q; then
// R = A - Q*B, and only its low deg-B coefficients matter. that is the whole
// divide: an inverse and two convolutions, O(deg A log deg A), no O(deg^2)
// schoolbook.
inline std::vector<i64> poly_mod(std::vector<i64> A, const std::vector<i64>& B) {
    trim(A);
    const int n = static_cast<int>(A.size()) - 1;  // deg A (-1 if A == 0).
    const int m = static_cast<int>(B.size()) - 1;  // deg B >= 1, B monic.
    if (n < m) return A;                            // already reduced.

    const int dq = n - m;  // deg Q.

    // reverse A (as degree n) and B (as degree m). B^R[0] = b_top = 1.
    std::vector<i64> Ar(A.rbegin(), A.rend());
    std::vector<i64> Br(B.rbegin(), B.rend());
    Ar.resize(static_cast<std::size_t>(dq + 1));  // only x^{dq+1} of A^R matters.

    // Q^R = A^R * (B^R)^{-1} mod x^{dq+1}, then flip back to Q (degree dq).
    const std::vector<i64> Bri = inverse(Br, dq + 1);
    std::vector<i64> Qr = convolve(Ar, Bri);
    Qr.resize(static_cast<std::size_t>(dq + 1), 0);
    const std::vector<i64> Q(Qr.rbegin(), Qr.rend());

    // R = A - Q*B, keep the low m coefficients (deg R < m).
    const std::vector<i64> QB = convolve(Q, B);
    std::vector<i64> R(static_cast<std::size_t>(m), 0);
    for (int i = 0; i < m; ++i) {
        const i64 ai = A[static_cast<std::size_t>(i)];
        const i64 qb = i < static_cast<int>(QB.size())
                           ? QB[static_cast<std::size_t>(i)]
                           : 0;
        i64 r = ai - qb;
        if (r < 0) r += kMod;  // ai, qb in [0, p) -> r in (-p, p), one add.
        R[static_cast<std::size_t>(i)] = r;
    }
    trim(R);
    return R;
}

// MULTIPOINT EVALUATION via the product / remainder tree.
//
// the naive road is Horner per point: O(n) each, O(n*m) total -- 1e10 at
// n=m=1e5, dead. the tree trades that for O(n log n + m log^2 m).
//
// LOWER BOUND. reading A is Omega(n), reading the points is Omega(m), and the
// output is m values so writing it is Omega(m). the product tree has O(log m)
// levels, each doing O(m log m) transform work across its nodes -> O(m log^2 m),
// plus one O(n log n) reduction of A against the root. you cannot beat the
// Omega(n+m) floor, and this sits a couple of logs above it.
//
// THE IDEA, one identity: A(x_j) == A(x) mod (x - x_j). the remainder of A on
// division by the linear factor (x - x_j) IS the value -- that's the remainder
// theorem. evaluating at all points is reducing A modulo every (x - x_j) at once,
// and the tree shares that work.
//
//   PRODUCT TREE (build). leaf j holds (x - x_j). each internal node holds the
//   product of its children -- one NTT multiply. the root holds prod_j (x - x_j).
//
//   REMAINDER TREE (descend). A mod (root) has degree < m and agrees with A at
//   every x_j (they differ by a multiple of prod (x - x_j), which is 0 at each
//   x_j). carry that remainder down: at a node with product P and children P_L,
//   P_R, pass R mod P_L left and R mod P_R right. reducing modulo a factor of P
//   preserves the value at every point that factor covers. at leaf j, R has
//   degree < 1 -- a constant, equal to A(x_j). done.
struct Evaluator {
    // segment-tree of node products, 1-indexed, sized to the classic 4*m bound.
    std::vector<std::vector<i64>> tree;
    std::vector<i64> pts;  // points, already reduced into [0, p).
    std::vector<i64> out;  // answers, filled leaf by leaf.

    // build node covering points [l, r). leaf = (x - x_l) = [-x_l, 1]; internal
    // = product of children. every node is monic -- product of monic linears --
    // which is exactly the invertibility poly_mod leans on.
    void build(int node, int l, int r) {
        if (r - l == 1) {
            const i64 x = pts[static_cast<std::size_t>(l)];
            tree[static_cast<std::size_t>(node)] = {x ? kMod - x : 0, 1};
            return;
        }
        const int mid = (l + r) >> 1;
        build(node << 1, l, mid);
        build((node << 1) | 1, mid, r);
        tree[static_cast<std::size_t>(node)] =
            convolve(tree[static_cast<std::size_t>(node << 1)],
                     tree[static_cast<std::size_t>((node << 1) | 1)]);
    }

    // descend with R = A mod (this node's product). at a leaf R is the constant
    // A(x_l); otherwise split R across the two child products.
    void eval(int node, int l, int r, std::vector<i64> R) {
        if (r - l == 1) {
            out[static_cast<std::size_t>(l)] = R.empty() ? 0 : R[0];
            return;
        }
        const int mid = (l + r) >> 1;
        const auto& L = tree[static_cast<std::size_t>(node << 1)];
        const auto& Rr = tree[static_cast<std::size_t>((node << 1) | 1)];
        eval(node << 1, l, mid, poly_mod(R, L));
        eval((node << 1) | 1, mid, r, poly_mod(std::move(R), Rr));
    }
};

// A(x_1), …, A(x_m). coefficients and points may arrive raw -- normalize into
// [0, p) up front, then the tree does the rest. m == 0 returns empty.
inline std::vector<i64> multipoint_eval(const std::vector<i64>& A_in,
                                        const std::vector<i64>& points) {
    const int m = static_cast<int>(points.size());
    if (m == 0) return {};

    Evaluator ev;
    ev.pts.resize(static_cast<std::size_t>(m));
    for (int j = 0; j < m; ++j) {
        i64 x = points[static_cast<std::size_t>(j)] % kMod;
        if (x < 0) x += kMod;
        ev.pts[static_cast<std::size_t>(j)] = x;
    }
    ev.out.assign(static_cast<std::size_t>(m), 0);
    ev.tree.assign(static_cast<std::size_t>(4 * m), {});
    ev.build(1, 0, m);

    // A reduced into [0, p). trailing zeros are fine -- poly_mod trims.
    std::vector<i64> A(A_in.size());
    for (std::size_t i = 0; i < A_in.size(); ++i) {
        i64 a = A_in[i] % kMod;
        if (a < 0) a += kMod;
        A[i] = a;
    }

    // one reduction against the root, then split down the tree.
    ev.eval(1, 0, m, poly_mod(std::move(A), ev.tree[1]));
    return ev.out;
}

}  // namespace p0052
