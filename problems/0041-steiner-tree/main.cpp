#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin:
//   n m k
//   m lines: u v w        (1-indexed vertices, undirected, w >= 1)
//   one line: k terminals (1-indexed)
// stdout: the minimum Steiner tree weight. 0 when k <= 1.
int main() {
    int n, m, k;
    if (std::scanf("%d %d %d", &n, &m, &k) != 3) return 1;

    std::vector<p0041::Edge> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) {
        int u, v;
        long long w;
        if (std::scanf("%d %d %lld", &u, &v, &w) != 3) return 1;
        edges.push_back({u - 1, v - 1, static_cast<std::int64_t>(w)});
    }

    std::vector<int> terminals(static_cast<std::size_t>(k));
    for (int i = 0; i < k; ++i) {
        int t;
        if (std::scanf("%d", &t) != 1) return 1;
        terminals[static_cast<std::size_t>(i)] = t - 1;
    }

    const std::int64_t ans = p0041::steiner_tree(n, edges, terminals);
    std::printf("%lld\n", static_cast<long long>(ans));
    return 0;
}
