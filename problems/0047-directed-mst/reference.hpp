#pragma once

#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "solution.hpp"  // reuse the Edge POD -- the algorithm below is its own.

namespace p0047 {

// dumb and obviously correct, for tiny n only. it does not contract anything. it
// just tries every way to hand each non-root vertex one of its incoming edges,
// throws out the ones that aren't a tree, and keeps the cheapest survivor.
//
//   an arborescence picks exactly one parent per non-root vertex. that's a
//   functional graph on the non-root vertices. with the root as the only
//   parentless node, it's a valid arborescence exactly when nobody loops -- every
//   parent chain has to terminate at the root, so the one thing to check is: no
//   cycle. sum the chosen weights, take the min over all cycle-free choices.
//
// cost is the product of the in-degrees, so this only survives for n <= 6..7 --
// which is the point. it's the ground truth the contraction is diffed against,
// including the cases that actually exercise the repair: forced cycles that must
// collapse, parallel edges where the min matters, and unreachable vertices that
// must come back -1.
//
// returns the minimum weight, or -1 if no arborescence exists.
inline std::int64_t min_arborescence_brute(int n,
                                           const std::vector<Edge>& edges,
                                           int root = 1) {
    if (n <= 0) return 0;
    const int r0 = root - 1;

    // incoming edges per vertex, as (parent, weight). the root takes none.
    std::vector<std::vector<std::pair<int, std::int64_t>>> incoming(
        static_cast<std::size_t>(n));
    for (const Edge& x : edges) {
        const int u = x.u - 1;
        const int v = x.v - 1;
        if (u == v) continue;    // self-loop -- a vertex can't be its own parent.
        if (u < 0 || u >= n || v < 0 || v >= n) continue;
        if (v == r0) continue;   // root is nobody's child.
        incoming[static_cast<std::size_t>(v)].emplace_back(u, x.w);
    }

    // a non-root vertex with nowhere to come from can never be covered -- no tree.
    for (int v = 0; v < n; ++v)
        if (v != r0 && incoming[static_cast<std::size_t>(v)].empty()) return -1;

    if (n == 1) return 0;

    std::vector<int> order;  // the non-root vertices, in a fixed choice order.
    for (int v = 0; v < n; ++v)
        if (v != r0) order.push_back(v);

    std::vector<int> parent(static_cast<std::size_t>(n), -1);
    std::int64_t best = std::numeric_limits<std::int64_t>::max();

    // does the chosen parent map reach the root from everywhere -- i.e. no cycle?
    auto acyclic = [&]() -> bool {
        for (const int start : order) {
            int x = start;
            int steps = 0;
            while (x != r0 && steps <= n) {
                x = parent[static_cast<std::size_t>(x)];
                ++steps;
            }
            if (x != r0) return false;  // ran past n hops -- it's looping.
        }
        return true;
    };

    // odometer over the in-edge lists: one choice per non-root vertex.
    std::function<void(std::size_t, std::int64_t)> rec =
        [&](std::size_t k, std::int64_t sum) {
            if (sum >= best) return;  // can't beat the incumbent -- prune.
            if (k == order.size()) {
                if (acyclic()) best = sum;
                return;
            }
            const int v = order[k];
            for (const auto& pw : incoming[static_cast<std::size_t>(v)]) {
                parent[static_cast<std::size_t>(v)] = pw.first;
                rec(k + 1, sum + pw.second);
            }
        };
    rec(0, 0);

    return best == std::numeric_limits<std::int64_t>::max() ? -1 : best;
}

}  // namespace p0047
