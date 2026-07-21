#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace p0047 {

// a directed edge u->v with weight w. weights are positive in the statement,
// but the contraction below reduces them, so the field is signed.
struct Edge {
    int u;
    int v;
    std::int64_t w;
};

// Minimum arborescence rooted at `root` (Chu-Liu/Edmonds).
//
// an arborescence is a spanning tree pointed away from the root: every non-root
// vertex owns exactly one incoming chosen edge, and following those edges
// backward from anywhere lands on the root. we want the cheapest one, or -1 when
// some vertex can't be reached from the root at all.
//
// the greedy instinct -- give each vertex its cheapest incoming edge -- is right
// until it isn't. those minimum in-edges can close a cycle, and a cycle is not a
// tree. that cycle is the whole problem, and it has no undirected analog:
// Prim/Kruskal grow one connected blob and never pick a *direction*, so a set of
// per-vertex minimum in-edges that loops back on itself is a shape they never
// produce and never have to repair. directed needs the repair. Chu-Liu is the
// repair.
//
// the repair is contraction. pick every vertex's minimum in-edge. no cycle among
// them -- done, sum is the answer. a cycle C -- pay for C's chosen in-edges up
// front, then collapse C into one supernode and ask the smaller question. the
// trick that makes the collapse honest: an edge e entering C at vertex v was
// competing against v's in-cycle edge, so charge it only the *extra* over that
// edge -- w(e) minus in[v]. we already paid in[v]. whatever cheaper edge the
// recursion later picks to enter the supernode, the books telescope to the exact
// tree weight. the reduced weight is w(e) - in[v] >= 0 because in[v] is the min
// over edges into v and e is one of them -- so nothing here ever goes negative.
//
// COMPLEXITY. each round is one linear sweep to pick min in-edges, one to find
// cycles, one to rebuild -- O(n + m). every round that finds a cycle merges at
// least two vertices into one, so the vertex count strictly drops and there are
// at most n rounds. O(n*m). reading the graph is already Omega(n + m), and each
// contraction removes >= 1 vertex so you get <= n of them -- that's the floor
// this rides.
//
// returns the minimum total weight, or -1 if the arborescence doesn't exist.
inline std::int64_t min_arborescence(int n, const std::vector<Edge>& edges_in,
                                     int root = 1) {
    if (n <= 0) return 0;
    const int r0 = root - 1;

    // 0-indexed working copy. self-loops carry no reachability and can only be
    // the minimum-in-edge trap for their own vertex, so drop them once here.
    std::vector<Edge> e;
    e.reserve(edges_in.size());
    for (const Edge& x : edges_in) {
        const int u = x.u - 1;
        const int v = x.v - 1;
        if (u == v) continue;                          // self-loop -- useless.
        if (u < 0 || u >= n || v < 0 || v >= n) continue;
        e.push_back(Edge{u, v, x.w});
    }

    // reachability first. an arborescence covers every vertex, so a vertex the
    // root can't reach kills the whole thing -- answer -1, before any work.
    {
        std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));
        for (const Edge& x : e)
            adj[static_cast<std::size_t>(x.u)].push_back(x.v);
        std::vector<char> seen(static_cast<std::size_t>(n), 0);
        std::vector<int> stk;
        stk.push_back(r0);
        seen[static_cast<std::size_t>(r0)] = 1;
        int reached = 1;
        while (!stk.empty()) {
            const int u = stk.back();
            stk.pop_back();
            for (const int w : adj[static_cast<std::size_t>(u)]) {
                if (!seen[static_cast<std::size_t>(w)]) {
                    seen[static_cast<std::size_t>(w)] = 1;
                    ++reached;
                    stk.push_back(w);
                }
            }
        }
        if (reached != n) return -1;
    }

    if (n == 1) return 0;  // lone root -- nothing to connect, weight zero.

    // reduced weights stay >= 0 and the answer is < n * max_w < 1e12, so int64
    // has all the room it needs. INF/4 keeps the "no in-edge" sentinel from
    // overflowing a comparison it never actually feeds into the sum.
    constexpr std::int64_t kInf = std::numeric_limits<std::int64_t>::max() / 4;

    std::int64_t res = 0;
    int cur_n = n;
    int cur_root = r0;

    while (true) {
        // cheapest incoming edge per vertex, and who it came from.
        std::vector<std::int64_t> in(static_cast<std::size_t>(cur_n), kInf);
        std::vector<int> pre(static_cast<std::size_t>(cur_n), -1);
        for (const Edge& x : e) {
            if (x.u != x.v && x.w < in[static_cast<std::size_t>(x.v)]) {
                in[static_cast<std::size_t>(x.v)] = x.w;
                pre[static_cast<std::size_t>(x.v)] = x.u;
            }
        }
        // a non-root supernode with no way in means an unreachable piece got
        // isolated by an earlier contraction -- no arborescence. the upfront BFS
        // catches this at round zero; this guards the contracted rounds.
        for (int i = 0; i < cur_n; ++i)
            if (i != cur_root && in[static_cast<std::size_t>(i)] == kInf)
                return -1;

        // walk the pre-pointers to label cycles. id != -1 marks a vertex already
        // placed in some cycle; vis == i marks a vertex seen on *this* walk, so
        // hitting it again is the cycle closing.
        std::vector<int> id(static_cast<std::size_t>(cur_n), -1);
        std::vector<int> vis(static_cast<std::size_t>(cur_n), -1);
        in[static_cast<std::size_t>(cur_root)] = 0;  // root pays nothing.
        int cnt = 0;
        for (int i = 0; i < cur_n; ++i) {
            res += in[static_cast<std::size_t>(i)];  // prepay every min in-edge.
            int v = i;
            while (vis[static_cast<std::size_t>(v)] != i &&
                   id[static_cast<std::size_t>(v)] == -1 && v != cur_root) {
                vis[static_cast<std::size_t>(v)] = i;
                v = pre[static_cast<std::size_t>(v)];
            }
            if (v != cur_root && id[static_cast<std::size_t>(v)] == -1) {
                for (int u = pre[static_cast<std::size_t>(v)]; u != v;
                     u = pre[static_cast<std::size_t>(u)])
                    id[static_cast<std::size_t>(u)] = cnt;
                id[static_cast<std::size_t>(v)] = cnt;
                ++cnt;
            }
        }

        if (cnt == 0) break;  // no cycle -- the prepaid sum is the answer.

        // vertices outside every cycle survive as their own singleton supernode.
        for (int i = 0; i < cur_n; ++i)
            if (id[static_cast<std::size_t>(i)] == -1)
                id[static_cast<std::size_t>(i)] = cnt++;

        // contract. an edge that crosses between supernodes gets charged only the
        // extra over the in-edge we already prepaid at its old head. edges inside
        // a supernode collapse to self-loops and drop out next round.
        for (Edge& x : e) {
            const int lv = x.v;
            x.u = id[static_cast<std::size_t>(x.u)];
            x.v = id[static_cast<std::size_t>(x.v)];
            if (x.u != x.v) x.w -= in[static_cast<std::size_t>(lv)];
        }
        cur_n = cnt;
        cur_root = id[static_cast<std::size_t>(cur_root)];
    }

    return res;
}

}  // namespace p0047
