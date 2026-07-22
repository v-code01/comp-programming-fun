#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 "n m"; then m lines "u v c w" -- an undirected edge u-v of color
//         c and weight w (vertices 1..n, colors 1..m, w any int; parallel edges,
//         multi-edges and self-loops allowed).
// stdout: one int64 -- the max total weight of an edge set that is a forest and
//         uses each color at most once. always >= 0 (the empty set is legal).
int main() {
    int n = 0, m = 0;
    if (std::scanf("%d %d", &n, &m) != 2) return 1;

    std::vector<std::tuple<int, int, int, std::int64_t>> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (int e = 0; e < m; ++e) {
        int u = 0, v = 0, c = 0;
        long long w = 0;
        if (std::scanf("%d %d %d %lld", &u, &v, &c, &w) != 4) return 1;
        edges.emplace_back(u, v, c, static_cast<std::int64_t>(w));
    }

    std::printf("%lld\n", static_cast<long long>(
                              p0053::max_weight_colorful_forest(n, edges)));
    return 0;
}
