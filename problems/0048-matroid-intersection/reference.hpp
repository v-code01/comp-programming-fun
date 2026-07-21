#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace p0048 {

// the dumb, obviously-correct oracle -- brute force over all 2^m edge subsets.
// no exchange graph, no matroid abstraction, not one line shared with the
// intersection. it just asks, for every subset of edges, two plain questions:
//   is it a forest?  (DSU: no edge closes a cycle)
//   is it colorful?  (no color repeats)
// keep the subsets that pass both, return the largest one's size. that IS the
// definition of a maximum common independent set, evaluated by exhaustion.
//
// O(2^m * m), so tiny m only -- this is the ground truth the augmenting-path
// solver is differential-tested against on thousands of small colored graphs.
inline int max_colorful_forest_ref(
    int n, const std::vector<std::tuple<int, int, int>>& edges) {
    const int m = static_cast<int>(edges.size());
    if (m == 0 || n <= 0) return 0;

    // pull endpoints and colors flat. 0-indexed vertices; find the palette size.
    std::vector<int> u(static_cast<std::size_t>(m)), v(static_cast<std::size_t>(m)),
        c(static_cast<std::size_t>(m));
    int max_color = 1;
    for (int j = 0; j < m; ++j) {
        const auto& [a, b, col] = edges[static_cast<std::size_t>(j)];
        u[static_cast<std::size_t>(j)] = a - 1;
        v[static_cast<std::size_t>(j)] = b - 1;
        c[static_cast<std::size_t>(j)] = col;
        if (col > max_color) max_color = col;
    }

    std::vector<int> parent(static_cast<std::size_t>(n));
    // stamped markers -- avoids clearing the color table once per subset.
    std::vector<unsigned> color_mark(static_cast<std::size_t>(max_color + 1), 0u);

    auto find = [&](int x) {
        while (parent[static_cast<std::size_t>(x)] != x)
            x = parent[static_cast<std::size_t>(x)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        return x;
    };

    int best = 0;
    const unsigned full = 1u << m;
    for (unsigned mask = 1; mask < full; ++mask) {
        const int cnt = __builtin_popcount(mask);
        if (cnt <= best) continue;  // can't beat the record -- don't even check.

        for (int i = 0; i < n; ++i) parent[static_cast<std::size_t>(i)] = i;
        bool ok = true;
        for (int j = 0; j < m && ok; ++j) {
            if (!((mask >> j) & 1u)) continue;
            // colorful: a color seen twice inside this subset kills it.
            const std::size_t col = static_cast<std::size_t>(c[static_cast<std::size_t>(j)]);
            if (color_mark[col] == mask + 1) {
                ok = false;
                break;
            }
            color_mark[col] = mask + 1;
            // forest: an edge whose endpoints already share a root closes a cycle.
            const int a = find(u[static_cast<std::size_t>(j)]);
            const int b = find(v[static_cast<std::size_t>(j)]);
            if (a == b) {
                ok = false;
                break;
            }
            parent[static_cast<std::size_t>(a)] = b;
        }
        if (ok) best = cnt;
    }
    return best;
}

}  // namespace p0048
