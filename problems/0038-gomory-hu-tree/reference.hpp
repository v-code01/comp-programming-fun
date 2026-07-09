#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

namespace p0038 {

// the dumb, obviously-correct oracle. no Gomory-Hu tree, no Dinic, no shared
// idea with the solution -- for one pair (u, v) it just runs a direct max-flow
// and returns it as the min cut (max-flow / min-cut). Edmonds-Karp: repeatedly
// BFS for the shortest augmenting path in the residual graph and push its
// bottleneck, until no path remains. residual lives in a dense n x n matrix, so
// undirected edges and summed parallels are a single symmetric cell.
//
// O(V * E^2) per pair -- glacial next to the tree, exact for tiny n. this is the
// ground truth the Gusfield build is diffed against, including disconnected
// graphs where the flow is 0 and the min cut has to read 0.
class DirectMaxFlow {
public:
    // build the dense symmetric capacity matrix once; every query reuses it by
    // copying it into a fresh residual. edges are 0-indexed {u, v, c}, parallels
    // summed into the same cell so the oracle sees the same graph the solver does.
    DirectMaxFlow(int n,
                  const std::vector<std::tuple<int, int, std::int64_t>>& edges)
        : n_(n),
          cap_(static_cast<std::size_t>(n),
               std::vector<std::int64_t>(static_cast<std::size_t>(n), 0)) {
        for (const auto& [u, v, c] : edges) {
            if (u == v) continue;  // self-loops never cross a cut.
            cap_[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)] += c;
            cap_[static_cast<std::size_t>(v)][static_cast<std::size_t>(u)] += c;
        }
    }

    std::int64_t min_cut(int s, int t) const {
        if (s == t) return 0;
        // fresh residual per query -- Edmonds-Karp mutates it down to zero flow.
        std::vector<std::vector<std::int64_t>> res = cap_;
        std::int64_t flow = 0;

        std::vector<int> parent(static_cast<std::size_t>(n_));
        std::vector<int> queue(static_cast<std::size_t>(n_));
        for (;;) {
            for (auto& p : parent) p = -1;
            parent[static_cast<std::size_t>(s)] = s;
            int qh = 0, qt = 0;
            queue[static_cast<std::size_t>(qt++)] = s;
            // BFS the residual for a shortest s->t augmenting path.
            while (qh < qt) {
                const int v = queue[static_cast<std::size_t>(qh++)];
                for (int u = 0; u < n_; ++u) {
                    if (parent[static_cast<std::size_t>(u)] == -1 &&
                        res[static_cast<std::size_t>(v)]
                           [static_cast<std::size_t>(u)] > 0) {
                        parent[static_cast<std::size_t>(u)] = v;
                        queue[static_cast<std::size_t>(qt++)] = u;
                    }
                }
            }
            if (parent[static_cast<std::size_t>(t)] == -1) break;  // t unreachable.

            // bottleneck along the path, then push it and credit the reverse.
            std::int64_t bottleneck = kInf;
            for (int v = t; v != s; v = parent[static_cast<std::size_t>(v)]) {
                const int p = parent[static_cast<std::size_t>(v)];
                const std::int64_t c =
                    res[static_cast<std::size_t>(p)][static_cast<std::size_t>(v)];
                if (c < bottleneck) bottleneck = c;
            }
            for (int v = t; v != s; v = parent[static_cast<std::size_t>(v)]) {
                const int p = parent[static_cast<std::size_t>(v)];
                res[static_cast<std::size_t>(p)][static_cast<std::size_t>(v)] -=
                    bottleneck;
                res[static_cast<std::size_t>(v)][static_cast<std::size_t>(p)] +=
                    bottleneck;
            }
            flow += bottleneck;
        }
        return flow;
    }

private:
    static constexpr std::int64_t kInf = 1'000'000'000'000'000'000LL;
    int n_;
    std::vector<std::vector<std::int64_t>> cap_;
};

}  // namespace p0038
