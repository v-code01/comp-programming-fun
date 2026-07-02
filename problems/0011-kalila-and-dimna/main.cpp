#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: n, then n heights a_i, then n recharge costs b_i.
// stdout: the minimum total recharge cost to fell every tree.
int main() {
    std::int64_t n;
    if (std::scanf("%lld", &n) != 1) return 1;

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    std::vector<std::int64_t> b(static_cast<std::size_t>(n));
    for (auto& x : a)
        if (std::scanf("%lld", &x) != 1) return 1;
    for (auto& x : b)
        if (std::scanf("%lld", &x) != 1) return 1;

    std::printf("%lld\n", static_cast<long long>(p0011::min_cost(a, b)));
    return 0;
}
