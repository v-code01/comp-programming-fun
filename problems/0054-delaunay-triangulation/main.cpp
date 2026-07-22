#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: n, then n lines "x y". stdout: sum of squared Delaunay edge lengths.
int main() {
    int n;
    if (std::scanf("%d", &n) != 1) return 1;

    std::vector<p0054::P> pts(static_cast<std::size_t>(n < 0 ? 0 : n));
    for (int i = 0; i < n; ++i) {
        long long x, y;
        if (std::scanf("%lld %lld", &x, &y) != 2) return 1;
        pts[static_cast<std::size_t>(i)] = {static_cast<std::int64_t>(x),
                                            static_cast<std::int64_t>(y)};
    }

    std::printf("%lld\n",
                static_cast<long long>(p0054::delaunay_edge_sq_sum(pts)));
    return 0;
}
