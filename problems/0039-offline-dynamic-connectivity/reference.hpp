#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

#include "solution.hpp"  // for p0039::Op -- shared 3-int op, nothing algorithmic.

namespace p0039 {

// the oracle. no segment tree, no rollback, no cleverness at all -- keep the ACTUAL
// current edge set, and for each query recompute connectivity from scratch. a
// std::set of normalized pairs is the whole state. adds insert, removes erase,
// queries flood-fill the live graph and check if v was reached from u.
//
// O(q * (n + q)) -- every query rebuilds adjacency and BFSes. slow on purpose. it
// shares nothing with the segtree solution but the answer it must agree with, which
// is exactly what makes it a trustworthy gate.
inline std::vector<int> referenceSolve(int n, const std::vector<Op>& ops) {
    std::set<std::pair<int, int>> edges;  // present edges, normalized (a < b).
    std::vector<int> answers;

    auto norm = [](int u, int v) {
        return u < v ? std::make_pair(u, v) : std::make_pair(v, u);
    };

    for (const Op& o : ops) {
        if (o.type == 1) {
            edges.insert(norm(o.u, o.v));
        } else if (o.type == 2) {
            edges.erase(norm(o.u, o.v));
        } else {
            // a vertex is always connected to itself -- no edge required.
            if (o.u == o.v) {
                answers.push_back(1);
                continue;
            }
            // build adjacency over the edges present right now, then BFS.
            std::vector<std::vector<int>> adj(static_cast<std::size_t>(n + 1));
            for (const auto& e : edges) {
                adj[static_cast<std::size_t>(e.first)].push_back(e.second);
                adj[static_cast<std::size_t>(e.second)].push_back(e.first);
            }
            std::vector<char> seen(static_cast<std::size_t>(n + 1), 0);
            std::vector<int> stack;
            stack.push_back(o.u);
            seen[static_cast<std::size_t>(o.u)] = 1;
            bool hit = false;
            while (!stack.empty() && !hit) {
                int x = stack.back();
                stack.pop_back();
                for (int w : adj[static_cast<std::size_t>(x)]) {
                    if (seen[static_cast<std::size_t>(w)]) continue;
                    if (w == o.v) {
                        hit = true;
                        break;
                    }
                    seen[static_cast<std::size_t>(w)] = 1;
                    stack.push_back(w);
                }
            }
            answers.push_back(hit ? 1 : 0);
        }
    }
    return answers;
}

}  // namespace p0039
