#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "solution.hpp"

// stdin: n / n costs / m / m coupon sizes. stdout: m lines, min money paid per
// coupon in input order.
int main() {
    int n;
    if (std::scanf("%d", &n) != 1) return 1;

    std::vector<std::int64_t> costs(static_cast<std::size_t>(n));
    for (auto& c : costs) {
        if (std::scanf("%lld", &c) != 1) return 1;
    }

    int m;
    if (std::scanf("%d", &m) != 1) return 1;

    // sort once up front, then each coupon is an O(1) lookup.
    p0002::Discount d(std::move(costs));
    for (int i = 0; i < m; ++i) {
        int q;
        if (std::scanf("%d", &q) != 1) return 1;
        std::printf("%lld\n", static_cast<long long>(d.answer(q)));
    }
    return 0;
}
