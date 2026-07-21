#include <cstdio>
#include <tuple>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 "n m"; then m lines "u v c" -- an undirected edge u-v of color c
//         (vertices 1..n, colors 1..m; parallel edges and multi-edges allowed).
// stdout: one int -- the max number of edges forming a forest that uses each
//         color at most once.
int main() {
    int n = 0, m = 0;
    if (std::scanf("%d %d", &n, &m) != 2) return 1;

    std::vector<std::tuple<int, int, int>> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (int e = 0; e < m; ++e) {
        int u = 0, v = 0, c = 0;
        if (std::scanf("%d %d %d", &u, &v, &c) != 3) return 1;
        edges.emplace_back(u, v, c);
    }

    std::printf("%d\n", p0048::max_colorful_forest(n, edges));
    return 0;
}
