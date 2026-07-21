#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

namespace p0048 {

// maximum colorful spanning forest -- MATROID INTERSECTION.
//
// two constraints on the same m edges, at once:
//   (A) forest -- no cycle. that's the GRAPHIC matroid M1.
//   (B) colorful -- at most one edge per color. that's a PARTITION matroid M2.
// pick the largest edge set independent in BOTH. the answer is |I|.
//
// why intersection is even tractable -- one matroid alone is a greedy problem,
// trivial. THREE matroids is NP-hard. two sits exactly on the boundary, and the
// reason is the exchange graph: one augmenting path per step, found by BFS, each
// step growing I by exactly one edge. that single path is the whole algorithm.
//
// the exchange graph D on the m edges, rebuilt every step around the current I:
//   sources  X1 = { y not in I : I+y is a forest }        -- free to add in M1.
//   sinks    X2 = { y not in I : color(y) unused in I }    -- free to add in M2.
//   arc x->y  (x in I, y not in I)  iff  I-x+y is a forest    -- M1 gates it.
//   arc y->x  (y not in I, x in I)  iff  I-x+y is colorful    -- M2 gates it.
// a shortest source-to-sink path P (fewest arcs, BFS) is the augmentation:
// I becomes I XOR V(P). the path alternates out,in,out,...,out -- one more
// out-node than in-node -- so it drops every in-node it touches and adds every
// out-node, and |I| climbs by exactly one. no path left means I is maximum.
//
// that's lawler's theorem: a MINIMUM-arc path keeps the symmetric difference
// independent in both matroids. an arbitrary augmenting path does NOT -- shortest
// is load-bearing. this is the trap that separates matroid intersection from
// plain bipartite matching, where any augmenting path works.
//
// the two gates, made cheap:
//   x->y (M1): if I+y is already a forest, I-x+y is a subset of it, independent
//     for EVERY x. but an arc INTO a source is dead weight -- the source already
//     sits at dist 0 -- so skip it. otherwise I+y closes exactly one cycle, its
//     fundamental cycle = the tree path between y's endpoints in the forest I.
//     removing x reopens it iff x lies ON that path. so x->y fires for exactly
//     the edges on that path.
//   y->x (M2): color(y) unused -> y is a sink, and a shortest path stops at the
//     first sink, so a sink's out-arcs never fire -- skip them. otherwise exactly
//     one edge in I holds color(y); the arc y->x fires only for x = that holder.
// dropping those two skippable arc classes can't shorten any real source->sink
// path, so BFS still returns a true shortest augmentation.
//
// COMPLEXITY -- one step is O(m*n): a DSU rebuild plus one forest-path scan per
// edge not in I. |I| <= n-1 (forest) and |I| <= #colors <= m (colorful), so the
// number of augmentations is <= r = min(n-1, m) -- each adds exactly one edge.
// total O(r*m*n), <= ~2.7e7 at n=m=300.
//
// LOWER BOUND -- reading the graph is Omega(n+m); you can't answer without
// touching every edge. above that, matroid intersection is polynomial -- the
// augmenting-path count is <= r and each path is polynomial to find. three
// matroids would be NP-hard; two is the tractable frontier, and the exchange
// graph is why.
class MatroidIntersection {
public:
    using Edge = std::tuple<int, int, int>;  // (u, v, color), u,v in 1..n.

    // 0 if there is nothing to pick. otherwise the size of a maximum common
    // independent set -- the answer.
    int solve(int n, const std::vector<Edge>& edges) {
        n_ = n;
        m_ = static_cast<int>(edges.size());
        if (m_ == 0 || n_ <= 0) return 0;

        eu_.resize(static_cast<std::size_t>(m_));
        ev_.resize(static_cast<std::size_t>(m_));
        ec_.resize(static_cast<std::size_t>(m_));
        int max_color = 1;
        for (int j = 0; j < m_; ++j) {
            const auto& [u, v, c] = edges[static_cast<std::size_t>(j)];
            eu_[static_cast<std::size_t>(j)] = u - 1;  // internal 0-indexed.
            ev_[static_cast<std::size_t>(j)] = v - 1;
            ec_[static_cast<std::size_t>(j)] = c;
            if (c > max_color) max_color = c;
        }
        max_color_ = max_color;

        in_i_.assign(static_cast<std::size_t>(m_), 0);
        while (augment()) {
            // each true is one edge more -- I grew by exactly one. false means no
            // source-to-sink path remains and I is a maximum common independent set.
        }

        int size = 0;
        for (int j = 0; j < m_; ++j) size += in_i_[static_cast<std::size_t>(j)];
        return size;
    }

private:
    int n_ = 0, m_ = 0, max_color_ = 1;
    std::vector<int> eu_, ev_, ec_;  // edge endpoints (0-indexed) and colors.
    std::vector<int> in_i_;          // membership of each edge in I.
    std::vector<int> dsu_;           // rebuilt every step over the edges in I.

    // path-halving find. the forest is tiny; correctness over cleverness.
    int find(int x) {
        while (dsu_[static_cast<std::size_t>(x)] != x) {
            dsu_[static_cast<std::size_t>(x)] =
                dsu_[static_cast<std::size_t>(dsu_[static_cast<std::size_t>(x)])];
            x = dsu_[static_cast<std::size_t>(x)];
        }
        return x;
    }

    // one augmentation: rebuild the exchange graph around I, BFS the shortest
    // source-to-sink path, XOR it into I. true iff I grew.
    bool augment() {
        // --- DSU over the current forest I. tells us fundamental cycles. ---
        dsu_.assign(static_cast<std::size_t>(n_), 0);
        for (int i = 0; i < n_; ++i) dsu_[static_cast<std::size_t>(i)] = i;
        // forest adjacency (only edges in I) for the tree-path scans below.
        std::vector<std::vector<std::pair<int, int>>> adj(
            static_cast<std::size_t>(n_));  // node -> (neighbor, edge index).
        // partition bookkeeping: which edge, if any, holds each color in I.
        std::vector<int> color_holder(static_cast<std::size_t>(max_color_ + 1), -1);
        for (int j = 0; j < m_; ++j) {
            if (!in_i_[static_cast<std::size_t>(j)]) continue;
            const int u = eu_[static_cast<std::size_t>(j)];
            const int v = ev_[static_cast<std::size_t>(j)];
            const int a = find(u), b = find(v);
            if (a != b) dsu_[static_cast<std::size_t>(a)] = b;  // I is a forest.
            adj[static_cast<std::size_t>(u)].push_back({v, j});
            adj[static_cast<std::size_t>(v)].push_back({u, j});
            color_holder[static_cast<std::size_t>(ec_[static_cast<std::size_t>(j)])] = j;
        }

        // --- arcs of the exchange graph, keyed by their tail node. ---
        std::vector<std::vector<int>> arc(static_cast<std::size_t>(m_));

        // reusable BFS buffers for the forest tree-path scans (stamped, no clear).
        std::vector<int> vis_stamp(static_cast<std::size_t>(n_), -1);
        std::vector<int> parent_v(static_cast<std::size_t>(n_), -1);
        std::vector<int> parent_e(static_cast<std::size_t>(n_), -1);
        int stamp = 0;
        std::vector<int> tree_q;

        for (int y = 0; y < m_; ++y) {
            if (in_i_[static_cast<std::size_t>(y)]) continue;  // out-nodes only.
            const int uy = eu_[static_cast<std::size_t>(y)];
            const int vy = ev_[static_cast<std::size_t>(y)];
            const bool is_source = (find(uy) != find(vy));
            const int cy = ec_[static_cast<std::size_t>(y)];
            const bool is_sink = (color_holder[static_cast<std::size_t>(cy)] == -1);

            // out->in (M2): a non-sink out-node reaches the sole holder of its color.
            // sink out-nodes are terminal in a shortest path, so skip their out-arcs.
            if (!is_sink)
                arc[static_cast<std::size_t>(y)].push_back(
                    color_holder[static_cast<std::size_t>(cy)]);

            // in->out (M1): the fundamental-cycle edges of a non-source out-node
            // point to it. source out-nodes sit at dist 0, so an arc into one can
            // never shorten a path -- skip building those.
            if (!is_source) {
                ++stamp;
                vis_stamp[static_cast<std::size_t>(uy)] = stamp;
                parent_v[static_cast<std::size_t>(uy)] = -1;
                parent_e[static_cast<std::size_t>(uy)] = -1;
                tree_q.clear();
                tree_q.push_back(uy);
                for (std::size_t h = 0; h < tree_q.size(); ++h) {
                    const int x = tree_q[h];
                    if (x == vy) break;
                    for (const auto& [nb, ei] : adj[static_cast<std::size_t>(x)]) {
                        if (vis_stamp[static_cast<std::size_t>(nb)] == stamp) continue;
                        vis_stamp[static_cast<std::size_t>(nb)] = stamp;
                        parent_v[static_cast<std::size_t>(nb)] = x;
                        parent_e[static_cast<std::size_t>(nb)] = ei;
                        tree_q.push_back(nb);
                    }
                }
                // walk the tree path uy..vy back, arc each edge x on it to y.
                for (int cur = vy; cur != uy;
                     cur = parent_v[static_cast<std::size_t>(cur)]) {
                    const int x = parent_e[static_cast<std::size_t>(cur)];
                    arc[static_cast<std::size_t>(x)].push_back(y);
                }
            }
        }

        // --- BFS: nearest sink from any source, over exchange-graph arcs. ---
        constexpr int kInf = 1 << 30;
        std::vector<int> dist(static_cast<std::size_t>(m_), kInf);
        std::vector<int> par(static_cast<std::size_t>(m_), -1);
        std::vector<int> bq;
        for (int y = 0; y < m_; ++y) {
            if (in_i_[static_cast<std::size_t>(y)]) continue;
            if (find(eu_[static_cast<std::size_t>(y)]) !=
                find(ev_[static_cast<std::size_t>(y)])) {
                dist[static_cast<std::size_t>(y)] = 0;
                bq.push_back(y);  // a source.
            }
        }

        int target = -1;
        for (std::size_t i = 0; i < bq.size() && target < 0; ++i) {
            const int c = bq[i];
            // out-node whose color is unused in I -> a sink. shortest path ends here.
            if (color_holder[static_cast<std::size_t>(ec_[static_cast<std::size_t>(c)])] ==
                    -1 &&
                !in_i_[static_cast<std::size_t>(c)]) {
                target = c;
                break;
            }
            for (const int nx : arc[static_cast<std::size_t>(c)]) {
                if (dist[static_cast<std::size_t>(nx)] != kInf) continue;
                dist[static_cast<std::size_t>(nx)] = dist[static_cast<std::size_t>(c)] + 1;
                par[static_cast<std::size_t>(nx)] = c;
                bq.push_back(nx);
            }
        }

        if (target < 0) return false;  // no augmenting path -- I is optimal.
        // XOR V(P) into I: flip membership along the whole path. adds the out-nodes,
        // drops the in-nodes -- one net edge gained.
        for (int cur = target; cur != -1; cur = par[static_cast<std::size_t>(cur)])
            in_i_[static_cast<std::size_t>(cur)] ^= 1;
        return true;
    }
};

// one-shot entry. edges are (u, v, color) with u,v in 1..n and color >= 1;
// parallel edges, multi-edges and self-loops are all accepted (a self-loop is
// simply never independent in M1, so it never enters I). returns the maximum
// number of edges that form a forest AND use each color at most once.
inline int max_colorful_forest(
    int n, const std::vector<std::tuple<int, int, int>>& edges) {
    MatroidIntersection mi;
    return mi.solve(n, edges);
}

}  // namespace p0048
