#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace p0014 {

// the oracle. dumb and obviously correct, small inputs only. keep a red/blue flag
// per vertex; a query BFSes out from v over the original tree and returns the
// distance to the first red vertex it reaches. O(n) per query, no cleverness --
// that's the point. the centroid solution has to agree with this everywhere.
//
// same API as NearestRed on purpose (build / paint / query) so a test can replay
// one op program on both and diff the answer streams.
class Bfs {
public:
    // vertex 0 starts red, matching the statement's initial condition -- both
    // solutions must fold it in identically or the differential is meaningless.
    void build(int n, const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        adj_.assign(static_cast<std::size_t>(n), {});
        for (const auto& [u, v] : edges) {
            adj_[static_cast<std::size_t>(u)].push_back(v);
            adj_[static_cast<std::size_t>(v)].push_back(u);
        }
        red_.assign(static_cast<std::size_t>(n), 0);
        red_[0] = 1;
    }

    void paint(int v) { red_[static_cast<std::size_t>(v)] = 1; }

    // BFS from v; the first red vertex popped is the nearest one. on a connected
    // tree with vertex 0 always red this always returns; -1 is unreachable.
    int query(int v) const {
        std::vector<int> dist(static_cast<std::size_t>(n_), -1);
        std::vector<int> q;
        q.push_back(v);
        dist[static_cast<std::size_t>(v)] = 0;
        for (std::size_t i = 0; i < q.size(); ++i) {
            const int u = q[i];
            if (red_[static_cast<std::size_t>(u)]) return dist[static_cast<std::size_t>(u)];
            for (const int w : adj_[static_cast<std::size_t>(u)]) {
                if (dist[static_cast<std::size_t>(w)] < 0) {
                    dist[static_cast<std::size_t>(w)] = dist[static_cast<std::size_t>(u)] + 1;
                    q.push_back(w);
                }
            }
        }
        return -1;
    }

private:
    int n_ = 0;
    std::vector<std::vector<int>> adj_;
    std::vector<char> red_;
};

}  // namespace p0014
