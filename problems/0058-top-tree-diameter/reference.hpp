#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0058 {

// the oracle. a weighted adjacency list and nothing clever. link and cut edit the
// lists; a diameter query floods the component and runs the textbook tree-diameter
// DP. O(component) per query, obviously correct, no clusters anywhere -- that's the
// point. the top tree has to agree with this on every op-sequence. same API
// (link/cut/diameter) so one program replays on both and the answers get diffed.
class Reference {
public:
    // vertices are 1..n; slot 0 is unused so indices match the top tree's 1-based ids.
    explicit Reference(int n) : n_(n), adj_(static_cast<std::size_t>(n + 1)) {}

    void link(int u, int v, std::int64_t w) {
        adj_[static_cast<std::size_t>(u)].push_back({v, w});
        adj_[static_cast<std::size_t>(v)].push_back({u, w});
    }

    // drop the edge from both endpoints' lists. it exists by contract.
    void cut(int u, int v) {
        erase(u, v);
        erase(v, u);
    }

    // diameter of u's tree: BFS the component, then a DFS DP. for each vertex keep the
    // two deepest downward reaches into distinct child subtrees; the diameter is the
    // best long1 + long2 over all vertices; a lone vertex is 0.
    std::int64_t diameter(int u) {
        std::vector<int> comp = component(u);
        // root the component at u and DFS iteratively (recursion could blow a 1e5 chain).
        // down[x] = deepest weighted reach from x into its subtree; ans folds in the two
        // deepest child reaches meeting at each vertex.
        std::vector<int> par(static_cast<std::size_t>(n_ + 1), 0);
        std::vector<int> order;
        order.reserve(comp.size());
        std::vector<char> seen(static_cast<std::size_t>(n_ + 1), 0);
        std::vector<int> stk{u};
        seen[static_cast<std::size_t>(u)] = 1;
        while (!stk.empty()) {
            int x = stk.back();
            stk.pop_back();
            order.push_back(x);
            for (auto [w, wt] : adj_[static_cast<std::size_t>(x)]) {
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                par[static_cast<std::size_t>(w)] = x;
                (void)wt;
                stk.push_back(w);
            }
        }
        std::vector<std::int64_t> down(static_cast<std::size_t>(n_ + 1), 0);
        std::int64_t ans = 0;
        // process children before parents -- reverse discovery order is a valid post-order.
        for (std::size_t i = order.size(); i-- > 0;) {
            int x = order[i];
            std::int64_t best1 = 0, best2 = 0;  // two deepest child reaches through x.
            for (auto [w, wt] : adj_[static_cast<std::size_t>(x)]) {
                // the forest is simple, so the parent is the one neighbor to skip.
                if (w == par[static_cast<std::size_t>(x)] && x != u) continue;
                std::int64_t reach = down[static_cast<std::size_t>(w)] + wt;
                if (reach > best1) {
                    best2 = best1;
                    best1 = reach;
                } else if (reach > best2) {
                    best2 = reach;
                }
            }
            down[static_cast<std::size_t>(x)] = best1;
            ans = std::max(ans, best1 + best2);  // longest path bending at x.
        }
        return ans;
    }

private:
    int n_;
    std::vector<std::vector<std::pair<int, std::int64_t>>> adj_;

    void erase(int a, int b) {
        auto& L = adj_[static_cast<std::size_t>(a)];
        for (std::size_t i = 0; i < L.size(); ++i)
            if (L[i].first == b) {
                L.erase(L.begin() + static_cast<std::ptrdiff_t>(i));
                return;
            }
    }

    std::vector<int> component(int u) {
        std::vector<char> seen(static_cast<std::size_t>(n_ + 1), 0);
        std::vector<int> q{u};
        seen[static_cast<std::size_t>(u)] = 1;
        for (std::size_t i = 0; i < q.size(); ++i)
            for (auto [w, wt] : adj_[static_cast<std::size_t>(q[i])]) {
                (void)wt;
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                q.push_back(w);
            }
        return q;
    }
};

}  // namespace p0058
