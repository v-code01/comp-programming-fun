#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0043 {

using i64 = std::int64_t;
using u64 = std::uint64_t;

// GF(p) with p = 998244353 = 119 * 2^23 + 1. that 2^23 is the whole reason this
// prime shows up everywhere -- p-1 is divisible by 2^23, so a primitive 2^k-th
// root of unity exists for every k up to 23. g = 3 generates the multiplicative
// group, so g^((p-1)/2^k) is a primitive 2^k-th root. that's the exact-arithmetic
// stand-in for e^{-2pi i / N} the complex FFT uses -- same butterflies, no doubles,
// no rounding, no precision argument to sweat.
inline constexpr i64 kMod = 998244353;
inline constexpr i64 kG = 3;  // primitive root of GF(998244353).

// (a*b) mod p. both factors sit under p < 2^30, so the product tops near 2^60 --
// it fits a u64 with room, no __int128 needed. and % by the constexpr kMod lowers
// to a magic-multiply, not a hardware divide, which is the difference in the NTT
// butterfly where this runs tens of millions of times. callers keep operands in
// [0, p); the u64 cast is safe for that whole range.
inline i64 mul(i64 a, i64 b) {
    return static_cast<i64>(static_cast<u64>(a) * static_cast<u64>(b) %
                            static_cast<u64>(kMod));
}

// a^e mod p, square-and-multiply. used for roots-of-unity and for the inverse.
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

// fermat inverse: p prime, so a^(p-2) == a^{-1} (mod p). the field is why the
// whole method exists -- newton divides by a_0 to start and the transform divides
// by the size to invert, and division only lives here because every nonzero
// residue has an inverse. over a composite modulus none of this holds.
inline i64 inv(i64 a) { return pw(a, kMod - 2); }

// iterative in-place radix-2 NTT. n must be a power of two, inputs already in
// [0, p). invert=false is the forward transform, invert=true the inverse (same
// butterflies with the inverse root, then a divide by n).
//
// three moving parts, same as the complex FFT:
//   bit-reversal permutation -- reorder so the in-place butterflies read
//     contiguous pairs. derived by walking a carry down the index.
//   butterfly with powers of the len-th root -- for each block of size len,
//     w = g^((p-1)/len) is a primitive len-th root; combine the low half with
//     w^j times the high half. inverse uses w^{-1}.
//   inverse scaling -- the forward-then-inverse round multiplies by n, so the
//     inverse pass finishes by multiplying every coefficient by n^{-1} mod p.
inline void ntt(std::vector<i64>& a, bool invert) {
    const int n = static_cast<int>(a.size());

    // bit-reversal permutation. j tracks the reverse of i by carrying the flip
    // from the top bit down -- no per-element ctz, just an incrementing mirror.
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // stage over block sizes len = 2, 4, …, n. each block folds its two halves.
    for (int len = 2; len <= n; len <<= 1) {
        // primitive len-th root of unity. len divides 2^23 divides p-1, so the
        // exponent (p-1)/len is an integer and g^that has exact order len.
        i64 w = pw(kG, (kMod - 1) / len);
        if (invert) w = inv(w);  // inverse transform walks the roots backward.
        const int half = len >> 1;
        for (int i = 0; i < n; i += len) {
            i64 wn = 1;  // w^0, w^1, … across the half-block. exact, no drift.
            for (int j = 0; j < half; ++j) {
                const i64 u = a[i + j];
                const i64 v = mul(a[i + j + half], wn);
                i64 s = u + v;
                if (s >= kMod) s -= kMod;  // u,v < p so u+v < 2p -- one subtract.
                i64 d = u - v;
                if (d < 0) d += kMod;  // u-v in (-p, p) -- one add.
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

// linear convolution over GF(p): res[k] = sum_i a[i] * b[k-i] mod p. pad both to
// a power of two large enough to hold the full product, transform, multiply
// pointwise, invert. exact -- the NTT is integer arithmetic end to end.
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

// B(x) = A(x)^{-1} mod x^n: the unique degree-<n polynomial with A*B == 1 mod x^n.
// requires a_0 != 0 (else A has no inverse -- a_0*b_0 == 1 is unsolvable).
//
// NEWTON DOUBLING. hold B correct mod x^k, push it to mod x^{2k} in one step:
//
//   B' = B * (2 - A*B)   (mod x^{2k}).
//
// WHY it doubles -- one line. suppose A*B == 1 (mod x^k). set E = 1 - A*B, so E is
// divisible by x^k (E == 0 mod x^k). then B' = B*(2 - A*B) = B*(1 + E), and
//
//   A*B' = (A*B)(2 - A*B) = (1 - E)(1 + E) = 1 - E^2.
//
// E == 0 mod x^k forces E^2 == 0 mod x^{2k}, so A*B' == 1 (mod x^{2k}). precision
// k -> 2k, every step. seed it at k=1 with B = [a_0^{-1}]: A*B == a_0*a_0^{-1} == 1
// (mod x^1). ceil(log2 n) doublings clear n.
//
// COST / LOWER BOUND, stated honestly. reading A is Omega(n). the answer has n
// coefficients, so writing it is Omega(n). newton runs O(log n) doublings; step at
// precision k truncates A to 2k terms and does O(1) convolutions of size O(k), so
// O(k log k) work. the k's double, 1, 2, 4, …, n, so the costs form a geometric
// series dominated by its last term: total O(n log n). the textbook O(n^2)
// recurrence (b_i = -a_0^{-1} sum_{j>=1} a_j b_{i-j}) is exactly what this beats --
// at n=1e5 that's 1e10 multiply-adds against ~n log n ~ 2e6-scale transform work.
inline std::vector<i64> inverse(const std::vector<i64>& A, int n) {
    // b_0 = a_0^{-1}. normalize a_0 into [0, p) first -- caller may hand raw input.
    i64 a0 = A[0] % kMod;
    if (a0 < 0) a0 += kMod;
    std::vector<i64> b{inv(a0)};

    // len is the precision B currently holds: A*B == 1 (mod x^len). double until it
    // covers n. the last step may overshoot past n -- harmless, we truncate at the
    // end, and truncating A to `t` terms each step keeps the series geometric.
    for (int len = 1; len < n; len <<= 1) {
        const int t = len << 1;  // target precision this step: mod x^t.

        // A truncated to x^t. terms of A at degree >= t never touch B mod x^t.
        std::vector<i64> at(A.begin(),
                            A.begin() + (static_cast<int>(A.size()) < t
                                             ? static_cast<int>(A.size())
                                             : t));

        // C = A*B, then D = 2 - C, both mod x^t.
        std::vector<i64> c = convolve(at, b);
        c.resize(static_cast<std::size_t>(t), 0);
        for (auto& x : c) x = x ? kMod - x : 0;  // negate: D = -C …
        c[0] += 2;                                // … + 2. c[0] < p+2, fix below.
        if (c[0] >= kMod) c[0] -= kMod;

        // B' = B*D mod x^t.
        std::vector<i64> nb = convolve(b, c);
        nb.resize(static_cast<std::size_t>(t), 0);
        b = std::move(nb);
    }

    b.resize(static_cast<std::size_t>(n), 0);
    return b;
}

}  // namespace p0043
