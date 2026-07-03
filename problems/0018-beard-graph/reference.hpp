#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace p0018 {

// the oracle. dumb and obviously correct, small inputs only. keep a color flag
// per edge; a query walks the a->b path edge by edge -- climb the deeper endpoint
// to its parent, checking that edge, until both endpoints meet at the LCA --
// counting edges and remembering whether it ever stepped over a white one. O(n)
// per query, no cleverness. the HLD solution has to agree with this everywhere.
//
// same API as BeardHLD on purpose (build / paint_black / paint_white / query) so
// a test can replay one op program on both and diff the answer streams.
class BeardWalk {
public:
    // vertex 0 is the root; every edge starts black -- matching the statement's
    // initial condition. both solutions must fold that in identically or the
    // differential is meaningless.
    void build(int n, const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        std::vector<std::vector<std::pair<int, int>>> adj(
            static_cast<std::size_t>(n));  // (neighbor, edge index).
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto [u, v] = edges[i];
            adj[static_cast<std::size_t>(u)].emplace_back(v, static_cast<int>(i));
            adj[static_cast<std::size_t>(v)].emplace_back(u, static_cast<int>(i));
        }

        parent_.assign(static_cast<std::size_t>(n), -1);
        depth_.assign(static_cast<std::size_t>(n), 0);
        up_edge_.assign(static_cast<std::size_t>(n), -1);  // edge to the parent.
        white_.assign(edges.size(), 0);                    // all black.

        // BFS from the root fixes parent, depth, and the edge above each vertex.
        std::vector<int> q;
        q.reserve(static_cast<std::size_t>(n));
        std::vector<char> seen(static_cast<std::size_t>(n), 0);
        q.push_back(0);
        seen[0] = 1;
        for (std::size_t i = 0; i < q.size(); ++i) {
            const int v = q[i];
            for (const auto [w, e] : adj[static_cast<std::size_t>(v)]) {
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                parent_[static_cast<std::size_t>(w)] = v;
                depth_[static_cast<std::size_t>(w)] =
                    depth_[static_cast<std::size_t>(v)] + 1;
                up_edge_[static_cast<std::size_t>(w)] = e;
                q.push_back(w);
            }
        }
    }

    void paint_black(int i) { white_[static_cast<std::size_t>(i)] = 0; }
    void paint_white(int i) { white_[static_cast<std::size_t>(i)] = 1; }

    // walk a and b up to their LCA, one edge at a time. count edges; the moment
    // (well, after the walk) any of them was white, the answer is -1.
    int query(int a, int b) const {
        if (a == b) return 0;
        int edges = 0;
        bool bad = false;

        // lift the deeper one until the depths match.
        while (depth_[static_cast<std::size_t>(a)] >
               depth_[static_cast<std::size_t>(b)]) {
            step(a, edges, bad);
        }
        while (depth_[static_cast<std::size_t>(b)] >
               depth_[static_cast<std::size_t>(a)]) {
            step(b, edges, bad);
        }
        // then lift both together until they collide at the LCA.
        while (a != b) {
            step(a, edges, bad);
            step(b, edges, bad);
        }
        return bad ? -1 : edges;
    }

private:
    // cross the edge above v: count it, note its color, move v to its parent.
    void step(int& v, int& edges, bool& bad) const {
        const int e = up_edge_[static_cast<std::size_t>(v)];
        if (white_[static_cast<std::size_t>(e)]) bad = true;
        ++edges;
        v = parent_[static_cast<std::size_t>(v)];
    }

    int n_ = 0;
    std::vector<int> parent_;
    std::vector<int> depth_;
    std::vector<int> up_edge_;  // input edge joining v to its parent, -1 at root.
    std::vector<char> white_;   // per-edge color, 1 = white.
};

}  // namespace p0018
