#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace p0053 {

// the dumb, obviously-correct oracle -- brute force over all 2^m edge subsets.
// no exchange graph, no matroid abstraction, no profit, not one line shared with
// the intersection. for every subset it asks the two plain questions:
//   is it a forest?  (DSU: no edge closes a cycle)
//   is it colorful?  (no color repeats)
// and if both pass, it scores the subset by its total weight and keeps the MAX.
// the empty set scores 0 and is always in the running, so the answer is >= 0 --
// that is the whole definition of the max-weight common independent set,
// evaluated by exhaustion.
//
// UNLIKE the unweighted twin, there is no popcount prune: a heavier subset can
// have FEWER edges (a lone +100 edge beats a colorful pair of +40s), so every
// subset has to be scored. O(2^m * m), tiny m only -- the ground truth the
// augmenting-path solver is differential-tested against on thousands of small
// weighted colored graphs.
inline std::int64_t max_weight_colorful_forest_ref(
    int n, const std::vector<std::tuple<int, int, int, std::int64_t>>& edges) {
    const int m = static_cast<int>(edges.size());
    if (m == 0 || n <= 0) return 0;

    // pull endpoints, colors, weights flat. 0-indexed vertices; palette size.
    std::vector<int> u(static_cast<std::size_t>(m)), v(static_cast<std::size_t>(m)),
        c(static_cast<std::size_t>(m));
    std::vector<std::int64_t> w(static_cast<std::size_t>(m));
    int max_color = 1;
    for (int j = 0; j < m; ++j) {
        const auto& [a, b, col, wt] = edges[static_cast<std::size_t>(j)];
        u[static_cast<std::size_t>(j)] = a - 1;
        v[static_cast<std::size_t>(j)] = b - 1;
        c[static_cast<std::size_t>(j)] = col;
        w[static_cast<std::size_t>(j)] = wt;
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

    std::int64_t best = 0;  // the empty set.
    const unsigned full = 1u << m;
    for (unsigned mask = 1; mask < full; ++mask) {
        for (int i = 0; i < n; ++i) parent[static_cast<std::size_t>(i)] = i;
        bool ok = true;
        std::int64_t weight = 0;
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
            weight += w[static_cast<std::size_t>(j)];
        }
        if (ok && weight > best) best = weight;
    }
    return best;
}

}  // namespace p0053
