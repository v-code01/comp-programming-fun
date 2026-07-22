#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0052 {

using i64 = std::int64_t;

// same prime as the solver, shared value and no shared logic. no NTT, no product
// tree, no division ever appears here. this is the ground truth the fast path
// answers to.
inline constexpr i64 kModRef = 998244353;

// tiny mod helpers, named apart from the solver's mul/pw so both headers coexist
// in one translation unit with no ODR clash. __int128 product so the multiply is
// obviously right -- correctness over cleverness, that is the whole job of the
// oracle.
inline i64 ref_mul(i64 a, i64 b) {
    return static_cast<i64>(static_cast<__int128>(a) * b % kModRef);
}

// ORACLE -- HORNER. evaluate A(x) = a_0 + a_1 x + … + a_{n-1} x^{n-1} at one
// point by nesting the multiplies: A(x) = ((… (a_{n-1} x + a_{n-2}) x + …) x) + a_0.
// walk the coefficients from the top down, one multiply-add each. O(n) per point.
// straight off the definition, no transform, no tree, no remainder -- a second,
// independent road to the same number.
inline i64 horner_eval(const std::vector<i64>& a, i64 x) {
    i64 xr = x % kModRef;
    if (xr < 0) xr += kModRef;
    i64 acc = 0;
    for (std::size_t k = a.size(); k-- > 0;) {
        i64 c = a[k] % kModRef;
        if (c < 0) c += kModRef;
        acc = (ref_mul(acc, xr) + c) % kModRef;
    }
    return acc;
}

// the whole answer the dumb way: Horner at every point. O(n*m) total -- the very
// cost the product tree exists to kill, but at small n, m it is instant and it
// cannot be subtly wrong the way a tree of NTTs can.
inline std::vector<i64> horner_multipoint(const std::vector<i64>& a,
                                          const std::vector<i64>& points) {
    std::vector<i64> out(points.size());
    for (std::size_t j = 0; j < points.size(); ++j)
        out[j] = horner_eval(a, points[j]);
    return out;
}

}  // namespace p0052
