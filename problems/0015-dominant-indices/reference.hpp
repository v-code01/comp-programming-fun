#pragma once

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace p0015 {

// the dumb, obviously-correct oracle -- O(n^2), small n only. no long paths, no
// shared pool, no offset arithmetic to be clever and wrong in. for each vertex v
// it walks v's subtree, tallies a fresh depth histogram relative to v, then reads
// off the argmax with the smallest depth on ties. that's it.
//
// the subtree of v is exactly the descendants of v in the tree rooted at 1. we
// root once (parent[] via BFS) and then, from v, only ever step to a neighbor w
// with parent[w] == current -- that never walks back up toward the root, so the
// traversal stays inside v's subtree.
//
// inputs and output match the fast solver: n >= 1, 1-indexed edges, out[i] is the
// dominant index of vertex i+1.
inline std::vector<int> dominant_indices_ref(
    int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<int> ans(static_cast<std::size_t>(n < 0 ? 0 : n), 0);
    if (n <= 1) return ans;

    const std::size_t N = static_cast<std::size_t>(n);

    // plain adjacency lists -- clarity over speed, the oracle only ever sees tiny
    // trees in the differential test.
    std::vector<std::vector<int>> g(N);
    for (const auto& e : edges) {
        g[static_cast<std::size_t>(e.first - 1)].push_back(e.second - 1);
        g[static_cast<std::size_t>(e.second - 1)].push_back(e.first - 1);
    }

    // root at 0, capture parent[] so subtree walks know which way is "down".
    std::vector<int> parent(N, -1);
    std::vector<int> order;
    order.reserve(N);
    std::vector<char> seen(N, 0);
    order.push_back(0);
    seen[0] = 1;
    for (std::size_t qi = 0; qi < order.size(); ++qi) {
        const int u = order[qi];
        for (const int w : g[static_cast<std::size_t>(u)]) {
            if (!seen[static_cast<std::size_t>(w)]) {
                seen[static_cast<std::size_t>(w)] = 1;
                parent[static_cast<std::size_t>(w)] = u;
                order.push_back(w);
            }
        }
    }

    std::vector<int> depth(N, 0);
    std::vector<int> stack;
    stack.reserve(N);
    std::vector<int> hist(N, 0);
    for (int v = 0; v < n; ++v) {
        // BFS the subtree of v, bucketing vertices by depth below v.
        std::fill(hist.begin(), hist.end(), 0);
        stack.clear();
        stack.push_back(v);
        depth[static_cast<std::size_t>(v)] = 0;
        int maxd = 0;
        for (std::size_t si = 0; si < stack.size(); ++si) {
            const int x = stack[si];
            const int d = depth[static_cast<std::size_t>(x)];
            ++hist[static_cast<std::size_t>(d)];
            if (d > maxd) maxd = d;
            for (const int w : g[static_cast<std::size_t>(x)]) {
                if (parent[static_cast<std::size_t>(w)] == x) {  // a child, step down.
                    depth[static_cast<std::size_t>(w)] = d + 1;
                    stack.push_back(w);
                }
            }
        }
        int best_d = 0, best_c = hist[0];
        for (int d = 1; d <= maxd; ++d) {
            if (hist[static_cast<std::size_t>(d)] > best_c) {
                best_c = hist[static_cast<std::size_t>(d)];
                best_d = d;
            }
        }
        ans[static_cast<std::size_t>(v)] = best_d;
    }
    return ans;
}

}  // namespace p0015
