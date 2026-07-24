#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace p0057 {

// k shortest s->t WALKS in a directed graph with non-negative weights, by
// Eppstein. a walk may repeat vertices and edges, so with a reachable cycle
// there are infinitely many -- we want the k smallest lengths in order.
//
// ---- the sidetrack decomposition (why this works) ----
// run Dijkstra on the REVERSED graph from t. dist[v] is the shortest v->t
// distance. pick, for every reachable v != t, one out-edge (v,x,w) with
// w + dist[x] == dist[v] -- that's v's TREE edge, and the tree edges form a
// shortest-path tree pointing at t.
//
// give every edge e=(u,v,w) a sidetrack cost
//     delta(e) = w + dist[v] - dist[u]  >= 0
// (a tree edge has delta 0; delta >= 0 because dist[u] is a min over out-edges).
// now decompose any s->t walk. walk along it; whenever you take a NON-tree edge
// you "sidetrack", and between sidetracks you ride tree edges straight downhill
// to the next sidetrack's tail. the walk is fully determined by the ORDERED list
// of sidetracks it takes, and telescoping the tree pieces gives
//     length(walk) = dist[s] + sum of delta over its sidetracks.
// so the shortest walk is dist[s] (no sidetracks), and every other walk is
// dist[s] plus a bag of sidetracks whose tails chain through the tree. the map
// walk <-> ordered sidetrack list is a bijection, so counting/ordering walks is
// counting/ordering these lists -- multiplicity and ties included.
//
// ---- the path graph (how we order the sidetrack lists cheaply) ----
// H_out(v): a min-heap (by delta) of v's sidetracks -- its out-edges minus the
// one tree edge (t keeps all of its, it has no tree edge).
// H_T(v): a PERSISTENT min-heap of every sidetrack you can still reach from v
// going downhill to t. it is H_T(tree_next(v)) with H_out(v) melded in, so the
// whole tree path v->t is shared structurally -- one meld per edge, nodes reused
// across every descendant. we use a persistent leftist heap for the sharing
// (the meldable-heap route; Eppstein's original 2-3-heap H_out shaves the build
// to O(m + n log n), we pay one more log on the m term and nothing on the query).
//
// best-first over the path graph. a state is "the ordered sidetrack list built
// so far", its key is the walk length. from a state whose last sidetrack is the
// heap node p (a heavier substitute or the next hop):
//   (a) move to p's heap child c in the same heap  -> swap that sidetrack for a
//       heavier sibling; length changes by delta[c]-delta[p] >= 0.
//   (b) cross to the root of H_T(head(p))            -> append the cheapest
//       sidetrack reachable after arriving at head(p); length grows by its delta.
// every edge weight is >= 0, so popping a global min-heap yields walk lengths in
// non-decreasing order. answer 0 is dist[s]; then each pop is the next length.
// each pop spawns <= 3 successors, so k answers cost O(k log k).
//
// ---- complexity / lower bound ----
// reading the graph is Omega(m). reverse Dijkstra is O(m + n log n). the build
// melds one sidetrack per edge into a persistent heap of depth O(log n), so
// O(m log n) (Eppstein's H_out variant: O(m + n log n)). the k answers come out
// in O(k log k). total O(m log n + n log n + k log k). the naive "grow every
// walk-prefix in a priority queue" (see reference.hpp) revisits shared structure
// and is worst-case exponential in the frontier -- Eppstein makes each of the k
// answers O(log k) by sharing H_T along tree paths and searching the path graph
// instead of the walk tree. lengths fit int64 by assumption; pushes that would
// overflow can't be among the k answers, so we drop them.

struct Edge {
    int u;
    int v;
    std::int64_t w;
};

// persistent leftist heap over sidetracks. a node carries the sidetrack's delta
// (the heap key) and the HEAD vertex of that edge (where its cross-edge lands).
// merge copies only the nodes on the right spine it touches, so old versions --
// the parents' H_T -- stay intact and shared.
struct SidetrackHeap {
    std::vector<int> left{0};
    std::vector<int> right{0};
    std::vector<int> rank{0};        // leftist "s-value"; rank[null]=0.
    std::vector<int> head{0};        // head vertex of the sidetrack edge.
    std::vector<std::int64_t> key{0};  // delta of the sidetrack.

    void reserve(std::size_t n) {
        left.reserve(n);
        right.reserve(n);
        rank.reserve(n);
        head.reserve(n);
        key.reserve(n);
    }

    // fresh singleton. never shared until a merge copies it.
    int make(std::int64_t k, int h) {
        left.push_back(0);
        right.push_back(0);
        rank.push_back(1);
        head.push_back(h);
        key.push_back(k);
        return static_cast<int>(left.size()) - 1;
    }

    // persistent meld. recursion rides the right spine only, so depth is
    // O(log size) and the explicit stack a leftist heap would need never grows.
    int merge(int a, int b) {
        if (a == 0) return b;
        if (b == 0) return a;
        if (key[static_cast<std::size_t>(a)] > key[static_cast<std::size_t>(b)])
            std::swap(a, b);
        // copy a -- everything below it must stay valid for older versions.
        const std::size_t as = static_cast<std::size_t>(a);
        const int c = static_cast<int>(left.size());
        left.push_back(left[as]);
        right.push_back(0);
        rank.push_back(0);
        head.push_back(head[as]);
        key.push_back(key[as]);
        const int merged = merge(right[as], b);
        const std::size_t cs = static_cast<std::size_t>(c);
        right[cs] = merged;
        if (rank[static_cast<std::size_t>(left[cs])] <
            rank[static_cast<std::size_t>(right[cs])])
            std::swap(left[cs], right[cs]);
        rank[cs] = rank[static_cast<std::size_t>(right[cs])] + 1;
        return c;
    }
};

// k shortest s->t walk lengths, non-decreasing. slot i is -1 if fewer than i+1
// walks exist. vertices are 1..n; edges directed with weight >= 0.
inline std::vector<std::int64_t> k_shortest_walks(
    int n, int s, int t, int k, const std::vector<Edge>& edges) {
    std::vector<std::int64_t> ans(static_cast<std::size_t>(std::max(k, 0)), -1);
    if (n <= 0 || k <= 0) return ans;

    constexpr std::int64_t INF = std::numeric_limits<std::int64_t>::max() / 4;
    const std::size_t m = edges.size();

    // ---- forward CSR: out-edges of each vertex, for the tree and the sidetracks.
    std::vector<int> fdeg(static_cast<std::size_t>(n) + 2, 0);
    for (const auto& e : edges) ++fdeg[static_cast<std::size_t>(e.u) + 1];
    for (int i = 1; i <= n + 1; ++i)
        fdeg[static_cast<std::size_t>(i)] += fdeg[static_cast<std::size_t>(i) - 1];
    std::vector<int> fhead = fdeg;
    std::vector<int> fto(m);
    std::vector<std::int64_t> fw(m);
    std::vector<int> epos(m);  // forward slot of edge e, so reverse can point back.
    for (std::size_t e = 0; e < m; ++e) {
        const auto& ed = edges[e];
        const std::size_t p =
            static_cast<std::size_t>(fhead[static_cast<std::size_t>(ed.u)]++);
        fto[p] = ed.v;
        fw[p] = ed.w;
        epos[e] = static_cast<int>(p);
    }
    fhead = fdeg;  // fhead[v]..fhead[v+1] is now v's out-edge range in fto/fw.

    // ---- reverse CSR: in-edges, so Dijkstra can flow from t outward. each
    // reverse arc remembers the forward slot of its edge, so the settling
    // relaxation can name the exact out-edge to treat as the tree edge. ----
    std::vector<int> rdeg(static_cast<std::size_t>(n) + 2, 0);
    for (const auto& e : edges) ++rdeg[static_cast<std::size_t>(e.v) + 1];
    for (int i = 1; i <= n + 1; ++i)
        rdeg[static_cast<std::size_t>(i)] += rdeg[static_cast<std::size_t>(i) - 1];
    std::vector<int> rhead = rdeg;
    std::vector<int> rto(m);
    std::vector<std::int64_t> rw(m);
    std::vector<int> rpos(m);
    for (std::size_t e = 0; e < m; ++e) {
        const auto& ed = edges[e];
        const std::size_t p =
            static_cast<std::size_t>(rhead[static_cast<std::size_t>(ed.v)]++);
        rto[p] = ed.u;
        rw[p] = ed.w;
        rpos[p] = epos[e];
    }
    rhead = rdeg;

    // ---- reverse Dijkstra from t. dist[v] = shortest v->t distance, and the
    // relaxation that settles v names v's TREE edge: parent = the neighbour that
    // handed v its final distance. reading the tree off Dijkstra (not off a
    // post-hoc "which out-edge realizes dist" scan) is what keeps it ACYCLIC even
    // with zero-weight edges -- a scan would happily pick a 0-weight self-loop,
    // which satisfies w + dist[x] == dist[v] with x == v and makes v its own
    // parent. the settle order can never do that: a parent is always popped
    // before its child. ----
    std::vector<std::int64_t> dist(static_cast<std::size_t>(n) + 1, INF);
    std::vector<int> tnext(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> tedge(static_cast<std::size_t>(n) + 1, -1);  // fto/fw index
    dist[static_cast<std::size_t>(t)] = 0;
    using DState = std::pair<std::int64_t, int>;  // (dist, vertex)
    std::priority_queue<DState, std::vector<DState>, std::greater<>> dpq;
    dpq.push({0, t});
    while (!dpq.empty()) {
        const auto [d, v] = dpq.top();
        dpq.pop();
        if (d > dist[static_cast<std::size_t>(v)]) continue;  // stale.
        for (int p = rhead[static_cast<std::size_t>(v)];
             p < rhead[static_cast<std::size_t>(v) + 1]; ++p) {
            const int u = rto[static_cast<std::size_t>(p)];
            const std::int64_t nd = d + rw[static_cast<std::size_t>(p)];
            if (nd < dist[static_cast<std::size_t>(u)]) {
                dist[static_cast<std::size_t>(u)] = nd;
                tnext[static_cast<std::size_t>(u)] = v;
                tedge[static_cast<std::size_t>(u)] = rpos[static_cast<std::size_t>(p)];
                dpq.push({nd, u});
            }
        }
    }

    if (dist[static_cast<std::size_t>(s)] >= INF) return ans;  // no s->t walk.
    ans[0] = dist[static_cast<std::size_t>(s)];  // the shortest walk.

    // ---- children lists in the tree (v hangs under tree_next(v)), so we can
    // build H_T parent-first with a BFS out from t. ----
    std::vector<int> cdeg(static_cast<std::size_t>(n) + 2, 0);
    for (int v = 1; v <= n; ++v) {
        if (v == t || dist[static_cast<std::size_t>(v)] >= INF) continue;
        ++cdeg[static_cast<std::size_t>(tnext[static_cast<std::size_t>(v)]) + 1];
    }
    for (int i = 1; i <= n + 1; ++i)
        cdeg[static_cast<std::size_t>(i)] += cdeg[static_cast<std::size_t>(i) - 1];
    std::vector<int> chead = cdeg;
    std::vector<int> cadj(
        static_cast<std::size_t>(cdeg[static_cast<std::size_t>(n) + 1]));
    for (int v = 1; v <= n; ++v) {
        if (v == t || dist[static_cast<std::size_t>(v)] >= INF) continue;
        const int par = tnext[static_cast<std::size_t>(v)];
        cadj[static_cast<std::size_t>(chead[static_cast<std::size_t>(par)]++)] = v;
    }
    chead = cdeg;

    // ---- build H_T(v) for every reachable vertex, parent before child. each
    // vertex melds its own sidetracks into the inherited (shared) parent heap. ----
    SidetrackHeap pool;
    pool.reserve(m + 1);  // grows past this on deep stacks; indices stay stable.
    std::vector<int> ht(static_cast<std::size_t>(n) + 1, 0);  // 0 == empty heap.

    std::vector<int> bfs;
    bfs.reserve(static_cast<std::size_t>(n));
    bfs.push_back(t);
    for (std::size_t qi = 0; qi < bfs.size(); ++qi) {
        const int v = bfs[qi];
        const int base = (v == t) ? 0 : ht[static_cast<std::size_t>(
                                             tnext[static_cast<std::size_t>(v)])];
        int h = base;
        for (int p = fhead[static_cast<std::size_t>(v)];
             p < fhead[static_cast<std::size_t>(v) + 1]; ++p) {
            if (p == tedge[static_cast<std::size_t>(v)]) continue;  // tree edge.
            const int x = fto[static_cast<std::size_t>(p)];
            if (dist[static_cast<std::size_t>(x)] >= INF) continue;  // dead head.
            const std::int64_t delta = fw[static_cast<std::size_t>(p)] +
                                       dist[static_cast<std::size_t>(x)] -
                                       dist[static_cast<std::size_t>(v)];
            // delta >= 0 always -- dist[v] is a min over v's out-edges.
            h = pool.merge(h, pool.make(delta, x));
        }
        ht[static_cast<std::size_t>(v)] = h;
        for (int c = chead[static_cast<std::size_t>(v)];
             c < chead[static_cast<std::size_t>(v) + 1]; ++c)
            bfs.push_back(cadj[static_cast<std::size_t>(c)]);
    }

    // ---- best-first over the path graph. state = (walk length, heap node). ----
    using PState = std::pair<std::int64_t, int>;  // (length, node)
    std::priority_queue<PState, std::vector<PState>, std::greater<>> pq;
    const int root = ht[static_cast<std::size_t>(s)];
    // len + add, but never overflow: an answer that fits int64 is guaranteed, so
    // a sum that would wrap is past every real answer and safe to drop.
    auto push = [&](std::int64_t len, std::int64_t add, int node) {
        if (node == 0 || add < 0) return;
        if (len > std::numeric_limits<std::int64_t>::max() - add) return;
        pq.push({len + add, node});
    };
    if (root != 0)
        push(ans[0], pool.key[static_cast<std::size_t>(root)], root);

    int cnt = 1;  // ans[0] is already the shortest walk.
    while (cnt < k && !pq.empty()) {
        const auto [len, node] = pq.top();
        pq.pop();
        ans[static_cast<std::size_t>(cnt++)] = len;
        const std::size_t nd = static_cast<std::size_t>(node);
        const std::int64_t kp = pool.key[nd];
        // (a) heavier siblings inside the same heap.
        const int lc = pool.left[nd];
        const int rc = pool.right[nd];
        if (lc) push(len, pool.key[static_cast<std::size_t>(lc)] - kp, lc);
        if (rc) push(len, pool.key[static_cast<std::size_t>(rc)] - kp, rc);
        // (b) append the cheapest sidetrack reachable after this edge's head.
        const int cross = ht[static_cast<std::size_t>(pool.head[nd])];
        if (cross) push(len, pool.key[static_cast<std::size_t>(cross)], cross);
    }

    return ans;
}

}  // namespace p0057
