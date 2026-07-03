#pragma once

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace p0013 {

// dumb and obviously correct, small trees only. build adjacency, root at 1 by a
// bfs to fix parents, then for every vertex walk its whole subtree from scratch,
// tally colors in a map, and sum the ones sitting at the max count. O(n) per
// vertex, O(n^2) overall -- the exact thing sack refuses to pay, kept here as the
// ground truth to diff against.
class ReferenceSolver {
public:
    std::vector<std::int64_t> solve(int n, const std::vector<int>& color,
                                    const std::vector<std::pair<int, int>>& edges) {
        if (n <= 0) return {};

        std::vector<std::vector<int>> g(static_cast<std::size_t>(n) + 1);
        for (const auto& e : edges) {
            g[static_cast<std::size_t>(e.first)].push_back(e.second);
            g[static_cast<std::size_t>(e.second)].push_back(e.first);
        }

        // root at 1: bfs to record each vertex's parent so subtree walks go down.
        std::vector<int> parent(static_cast<std::size_t>(n) + 1, 0);
        std::vector<int> order;  // bfs order, also a valid top-down sweep.
        order.reserve(static_cast<std::size_t>(n));
        std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
        std::vector<int> queue;
        queue.push_back(1);
        seen[1] = 1;
        parent[1] = 0;
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            int v = queue[qi];
            order.push_back(v);
            for (int u : g[static_cast<std::size_t>(v)]) {
                if (seen[static_cast<std::size_t>(u)]) continue;
                seen[static_cast<std::size_t>(u)] = 1;
                parent[static_cast<std::size_t>(u)] = v;
                queue.push_back(u);
            }
        }

        std::vector<std::int64_t> answer(static_cast<std::size_t>(n) + 1, 0);
        for (int root = 1; root <= n; ++root) {
            // walk subtree(root) with an explicit stack, never stepping to parent.
            std::unordered_map<int, int> cnt;
            std::vector<int> stack;
            stack.push_back(root);
            int max_cnt = 0;
            std::int64_t sum_at_max = 0;
            while (!stack.empty()) {
                int v = stack.back();
                stack.pop_back();
                int col = color[static_cast<std::size_t>(v - 1)];
                int c = ++cnt[col];
                if (c > max_cnt) {
                    max_cnt = c;
                    sum_at_max = col;
                } else if (c == max_cnt) {
                    sum_at_max += col;
                }
                for (int u : g[static_cast<std::size_t>(v)]) {
                    if (u == parent[static_cast<std::size_t>(v)]) continue;
                    stack.push_back(u);
                }
            }
            answer[static_cast<std::size_t>(root)] = sum_at_max;
        }

        std::vector<std::int64_t> out(static_cast<std::size_t>(n));
        for (int v = 1; v <= n; ++v)
            out[static_cast<std::size_t>(v - 1)] = answer[static_cast<std::size_t>(v)];
        return out;
    }
};

inline std::vector<std::int64_t> reference_solve(
    int n, const std::vector<int>& color,
    const std::vector<std::pair<int, int>>& edges) {
    ReferenceSolver r;
    return r.solve(n, color, edges);
}

}  // namespace p0013
