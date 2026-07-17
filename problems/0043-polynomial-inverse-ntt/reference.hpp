#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0043 {

using i64 = std::int64_t;

// same prime as the solver -- shared value, shared no logic. the NTT and newton
// never appear here. this is the ground truth they answer to.
inline constexpr i64 kModRef = 998244353;

// tiny mod helpers, named apart from the solver's so both headers coexist in one
// translation unit without an ODR clash. obviously correct, that's the point.
inline i64 ref_mul(i64 a, i64 b) {
    return static_cast<i64>(static_cast<__int128>(a) * b % kModRef);
}
inline i64 ref_pw(i64 a, i64 e) {
    a %= kModRef;
    if (a < 0) a += kModRef;
    i64 r = 1;
    while (e) {
        if (e & 1) r = ref_mul(r, a);
        a = ref_mul(a, a);
        e >>= 1;
    }
    return r;
}
inline i64 ref_inv(i64 a) { return ref_pw(a, kModRef - 2); }  // fermat.

// ORACLE (a) -- the direct O(n^2) inverse, straight off the definition. match
// A(x)*B(x) == 1 (mod x^n) coefficient by coefficient. degree 0 gives a_0 b_0 = 1,
// so b_0 = a_0^{-1}. degree i >= 1 gives sum_{j=0..i} a_j b_{i-j} = 0, i.e.
//
//   a_0 b_i + sum_{j=1..i} a_j b_{i-j} = 0  ->  b_i = -a_0^{-1} sum_{j=1..i} a_j b_{i-j}.
//
// no NTT, no newton -- just the recurrence, O(n^2). requires a_0 != 0.
inline std::vector<i64> direct_inverse(const std::vector<i64>& A, int n) {
    std::vector<i64> a(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n && i < static_cast<int>(A.size()); ++i) {
        i64 v = A[static_cast<std::size_t>(i)] % kModRef;
        if (v < 0) v += kModRef;
        a[static_cast<std::size_t>(i)] = v;
    }

    std::vector<i64> b(static_cast<std::size_t>(n), 0);
    const i64 a0inv = ref_inv(a[0]);
    b[0] = a0inv;
    for (int i = 1; i < n; ++i) {
        i64 s = 0;
        for (int j = 1; j <= i; ++j)
            s = (s + ref_mul(a[static_cast<std::size_t>(j)],
                             b[static_cast<std::size_t>(i - j)])) %
                kModRef;
        // b_i = -a_0^{-1} * s.
        b[static_cast<std::size_t>(i)] =
            ref_mul(a0inv, (kModRef - s) % kModRef);
    }
    return b;
}

// ORACLE (b) -- the verifier. brute O(n^2) convolution of A and a candidate B,
// truncated to x^n, must equal [1, 0, 0, …, 0]. this checks the answer without
// ever computing it: whatever produced B, A*B mod x^n either is 1 or it isn't.
// independent of oracle (a) -- a second road to the same verdict.
inline bool verify_inverse(const std::vector<i64>& A, const std::vector<i64>& B,
                           int n) {
    for (int k = 0; k < n; ++k) {
        i64 s = 0;
        for (int j = 0; j <= k; ++j) {
            if (j >= static_cast<int>(A.size()) ||
                (k - j) >= static_cast<int>(B.size()))
                continue;
            i64 aj = A[static_cast<std::size_t>(j)] % kModRef;
            if (aj < 0) aj += kModRef;
            i64 bk = B[static_cast<std::size_t>(k - j)] % kModRef;
            if (bk < 0) bk += kModRef;
            s = (s + ref_mul(aj, bk)) % kModRef;
        }
        const i64 want = (k == 0) ? 1 : 0;
        if (s != want) return false;
    }
    return true;
}

}  // namespace p0043
