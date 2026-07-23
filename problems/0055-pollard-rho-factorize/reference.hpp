#pragma once

#include <cstdint>
#include <vector>

namespace p0055 {

// trial division -- dumb and obviously correct, the thing we trust. divide out 2,
// then every odd d up to sqrt(n). each divisor that hits is prime by
// construction, since everything smaller is already gone. whatever survives
// above sqrt(n) is itself a prime factor, so push it. the list comes out
// non-decreasing for free.
//
// O(sqrt n): exact up to ~1e12 (sqrt = 1e6), hopeless past that. no rho, no
// miller-rabin -- that's the point. it shares no code and no idea with the fast
// path, so agreement between them is real evidence, not a tautology.
inline std::vector<std::uint64_t> trial_factor(std::uint64_t n) {
    std::vector<std::uint64_t> f;
    while (n > 1 && (n & 1) == 0) {
        f.push_back(2);
        n >>= 1;
    }
    // d <= n / d instead of d * d <= n -- same test, no overflow of d*d near 1e12.
    for (std::uint64_t d = 3; d <= n / d; d += 2) {
        while (n % d == 0) {
            f.push_back(d);
            n /= d;
        }
    }
    if (n > 1) f.push_back(n);
    return f;
}

}  // namespace p0055
