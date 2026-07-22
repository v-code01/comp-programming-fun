#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: n, then n lines "x y z". stdout: number of triangular hull faces.
int main() {
    int n;
    if (std::scanf("%d", &n) != 1) return 1;

    std::vector<p0051::P> pts(static_cast<std::size_t>(n < 0 ? 0 : n));
    for (int i = 0; i < n; ++i) {
        long long x, y, z;
        if (std::scanf("%lld %lld %lld", &x, &y, &z) != 3) return 1;
        pts[static_cast<std::size_t>(i)] = {static_cast<std::int64_t>(x),
                                            static_cast<std::int64_t>(y),
                                            static_cast<std::int64_t>(z)};
    }

    std::printf("%d\n", p0051::convex_hull_faces(pts));
    return 0;
}
