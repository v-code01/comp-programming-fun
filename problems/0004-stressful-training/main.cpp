#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: "n k", then n starting charges a_i, then n drains b_i.
// stdout: the smallest charger power that saves everyone, or -1.
int main() {
    std::int64_t n, k;
    if (std::scanf("%lld %lld", &n, &k) != 2) return 1;

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    std::vector<std::int64_t> b(static_cast<std::size_t>(n));
    for (auto& x : a)
        if (std::scanf("%lld", &x) != 1) return 1;
    for (auto& x : b)
        if (std::scanf("%lld", &x) != 1) return 1;

    std::printf("%lld\n", static_cast<long long>(p0004::min_charger_power(a, b, k)));
    return 0;
}
