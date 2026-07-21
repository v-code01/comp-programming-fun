#pragma once

#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace p0046 {

// one op in the statement's alphabet, vertices 1-based.
//   type 1: insert edge u--v (absent by contract)
//   type 2: delete edge u--v (present by contract)
//   type 3: query connected(u, v) -> 1 / 0
struct Op {
    int type;
    int u;
    int v;
};

// the oracle. no forest, no levels, no treap -- keep the ACTUAL live edge set and
// recompute connectivity from scratch on every query. a std::set of normalized
// pairs is the whole state. inserts add, deletes erase, queries flood-fill the
// current graph and check whether v is reachable from u.
//
// O(q * (n + q)) -- every query rebuilds adjacency and BFSes. slow on purpose. it
// shares nothing with the HDLT structure but the answer stream it must match,
// which is exactly what makes it a trustworthy gate. unlike the offline oracle in
// 0039 it consumes ops one at a time, in order -- same online contract as the
// structure it checks.
inline std::vector<int> referenceSolve(int n, const std::vector<Op>& ops) {
    std::set<std::pair<int, int>> edges;  // live edges, normalized (a < b).
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
            if (o.u == o.v) {  // a vertex is always connected to itself.
                answers.push_back(1);
                continue;
            }
            std::vector<std::vector<int>> adj(static_cast<std::size_t>(n) + 1);
            for (const auto& e : edges) {
                adj[static_cast<std::size_t>(e.first)].push_back(e.second);
                adj[static_cast<std::size_t>(e.second)].push_back(e.first);
            }
            std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
            std::vector<int> stack;
            stack.push_back(o.u);
            seen[static_cast<std::size_t>(o.u)] = 1;
            bool hit = false;
            while (!stack.empty() && !hit) {
                int x = stack.back();
                stack.pop_back();
                for (int w : adj[static_cast<std::size_t>(x)]) {
                    if (w == o.v) {
                        hit = true;
                        break;
                    }
                    if (seen[static_cast<std::size_t>(w)]) continue;
                    seen[static_cast<std::size_t>(w)] = 1;
                    stack.push_back(w);
                }
            }
            answers.push_back(hit ? 1 : 0);
        }
    }
    return answers;
}

}  // namespace p0046
