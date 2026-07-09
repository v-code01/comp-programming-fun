#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

namespace p0038 {

// all-pairs minimum cut on an undirected weighted graph -- the Gomory-Hu tree.
//
// the min cut between u and v equals the max flow between u and v (max-flow /
// min-cut, Ford & Fulkerson 1956). the naive answer runs a max-flow for every
// one of the C(n,2) pairs. Gomory & Hu (1961) proved you never need more than
// n-1. Gusfield (1990) then stripped the construction down to something you can
// hold in your head: n-1 max-flows on the ORIGINAL graph, no vertex contraction.
//
// ---- lower bound, argued honestly ----
//
// reading the graph is Omega(n + m) -- you can't answer a cut query about an
// edge you never looked at. and the build needs Omega(n-1) max-flow calls: the
// tree has n-1 edges, each carrying a distinct cut value, and each value is the
// outcome of a genuine s-t max-flow. so n-1 max-flows isn't a heuristic budget,
// it's the count the tree structure forces. the win over naive is exact: n-1
// flows instead of n(n-1)/2, and the resulting tree encodes ALL pairwise cuts.
// that's the whole point of Gomory-Hu -- one small tree answers every pair.
//
// ---- why a tree can hold every pairwise cut ----
//
// Gomory-Hu: there's a weighted tree on the same n vertices where, for any u, v,
//   mincut(u, v) = the MINIMUM edge weight on the u-v path in the tree.
// intuition: the lightest edge on the path is a bottleneck you can't route
// around inside the tree, and the theorem says the graph can't route around it
// either. cut that one tree edge and the two sides are exactly a min u-v cut in
// the graph. so all-pairs min cut collapses to all-pairs path-minimum on a tree.
//
// disconnected graphs fall out for free: if u and v sit in different components
// their max flow is 0, that 0 lands on a tree edge, and the path-min through it
// reads 0 -- no special case in the query path.

// Dinic's max-flow, resettable so one built graph serves all n-1 runs. undirected
// edges: an edge {u, v, c} becomes two residual arcs BOTH starting at capacity c
// (forward and twin), so flow may push either way -- that's what makes this an
// undirected cut. O(V^2 * E) worst case, far quicker here; n <= 150 keeps every
// run cheap. capacities are int64: a vertex can carry ~149 * 1e9 ~ 1.5e11 of
// incident capacity, past 32 bits, and the running flow lives at that scale too.
class Dinic {
public:
    explicit Dinic(int num_nodes)
        : head_(static_cast<std::size_t>(num_nodes), -1),
          level_(static_cast<std::size_t>(num_nodes)),
          iter_(static_cast<std::size_t>(num_nodes)),
          n_(num_nodes) {}

    // one undirected edge = two arcs, each at full capacity c. twin of arc e is
    // e^1, so pushing on one and crediting the other is a single xor, no lookup.
    void add_undirected(int u, int v, std::int64_t c) {
        to_.push_back(v);
        cap_.push_back(c);
        next_.push_back(head_[static_cast<std::size_t>(u)]);
        head_[static_cast<std::size_t>(u)] = static_cast<int>(to_.size()) - 1;

        to_.push_back(u);
        cap_.push_back(c);  // undirected -- the reverse arc also starts at c.
        next_.push_back(head_[static_cast<std::size_t>(v)]);
        head_[static_cast<std::size_t>(v)] = static_cast<int>(to_.size()) - 1;

        orig_.push_back(c);  // remember both arcs' starting caps for the reset.
        orig_.push_back(c);
    }

    // wind every arc back to its start capacity. Gusfield runs a fresh max-flow
    // between i and parent[i] on the SAME original graph each round -- residual
    // left over from the last run would poison the next cut.
    void reset() { cap_ = orig_; }

    std::int64_t max_flow(int s, int t) {
        std::int64_t flow = 0;
        while (bfs(s, t)) {
            for (int v = 0; v < n_; ++v) {
                iter_[static_cast<std::size_t>(v)] =
                    head_[static_cast<std::size_t>(v)];
            }
            std::int64_t f;
            while ((f = dfs(s, t, kInf)) > 0) flow += f;
        }
        return flow;
    }

    // the min-cut side: vertices reachable from s in the residual graph after a
    // max-flow. this IS the S side of a minimum s-t cut -- the arcs leaving it
    // are saturated, and their capacity sums to the flow value. Gusfield reads
    // exactly this set to decide which subtree a vertex belongs to.
    std::vector<char> min_cut_side(int s) const {
        std::vector<char> seen(static_cast<std::size_t>(n_), 0);
        std::vector<int> stack;
        stack.push_back(s);
        seen[static_cast<std::size_t>(s)] = 1;
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            for (int e = head_[static_cast<std::size_t>(v)]; e != -1;
                 e = next_[static_cast<std::size_t>(e)]) {
                const int u = to_[static_cast<std::size_t>(e)];
                if (cap_[static_cast<std::size_t>(e)] > 0 &&
                    !seen[static_cast<std::size_t>(u)]) {
                    seen[static_cast<std::size_t>(u)] = 1;
                    stack.push_back(u);
                }
            }
        }
        return seen;
    }

private:
    static constexpr std::int64_t kInf = 1'000'000'000'000'000'000LL;

    bool bfs(int s, int t) {
        for (auto& l : level_) l = -1;
        level_[static_cast<std::size_t>(s)] = 0;
        bfs_q_.clear();
        bfs_q_.push_back(s);
        for (std::size_t qi = 0; qi < bfs_q_.size(); ++qi) {
            const int v = bfs_q_[qi];
            for (int e = head_[static_cast<std::size_t>(v)]; e != -1;
                 e = next_[static_cast<std::size_t>(e)]) {
                const int u = to_[static_cast<std::size_t>(e)];
                if (cap_[static_cast<std::size_t>(e)] > 0 &&
                    level_[static_cast<std::size_t>(u)] < 0) {
                    level_[static_cast<std::size_t>(u)] =
                        level_[static_cast<std::size_t>(v)] + 1;
                    bfs_q_.push_back(u);
                }
            }
        }
        return level_[static_cast<std::size_t>(t)] >= 0;
    }

    std::int64_t dfs(int v, int t, std::int64_t f) {
        if (v == t) return f;
        for (int& e = iter_[static_cast<std::size_t>(v)]; e != -1;
             e = next_[static_cast<std::size_t>(e)]) {
            const int u = to_[static_cast<std::size_t>(e)];
            if (cap_[static_cast<std::size_t>(e)] > 0 &&
                level_[static_cast<std::size_t>(u)] ==
                    level_[static_cast<std::size_t>(v)] + 1) {
                const std::int64_t d =
                    dfs(u, t,
                        f < cap_[static_cast<std::size_t>(e)]
                            ? f
                            : cap_[static_cast<std::size_t>(e)]);
                if (d > 0) {
                    cap_[static_cast<std::size_t>(e)] -= d;
                    cap_[static_cast<std::size_t>(e ^ 1)] += d;
                    return d;
                }
            }
        }
        return 0;
    }

    std::vector<int> head_;
    std::vector<int> to_;
    std::vector<std::int64_t> cap_;
    std::vector<std::int64_t> orig_;
    std::vector<int> next_;
    std::vector<int> level_;
    std::vector<int> iter_;
    std::vector<int> bfs_q_;
    int n_;
};

// build the Gomory-Hu tree with Gusfield's method, then answer any pair in O(1)
// from a precomputed path-minimum table. edges are 0-indexed {u, v, c}; parallel
// edges should be summed by the caller before they arrive here.
class GomoryHu {
public:
    GomoryHu(int n, const std::vector<std::tuple<int, int, std::int64_t>>& edges)
        : n_(n),
          parent_(static_cast<std::size_t>(n), 0),
          fl_(static_cast<std::size_t>(n), 0) {
        build(edges);
        precompute_pairs();
    }

    // min cut between u and v = lightest tree edge on their path. answered from
    // the table, so it's a single load.
    std::int64_t min_cut(int u, int v) const {
        return answer_[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)];
    }

private:
    // Gusfield's construction, tree rooted at 0. hold parent[i] = i's neighbour
    // toward the root and fl[i] = the cut on that tree edge. for each s in 1..n-1
    // cut s from its current parent t on the ORIGINAL graph; the flow value is
    // the tree-edge weight, and the min-cut side X (vertices on s's side) tells
    // us which other vertices should now hang off s instead of t.
    //
    // the grandparent rotation is the subtle half. if t's own parent landed on
    // s's side of the cut, then s sits closer to the root than t does -- so swap
    // them: s takes t's parent, t hangs under s, and the two edge weights trade.
    // skip that check and vertices can end up on the wrong subtree, which is the
    // bug the differential oracle exists to catch.
    void build(const std::vector<std::tuple<int, int, std::int64_t>>& edges) {
        Dinic din(n_);
        for (const auto& [u, v, c] : edges) din.add_undirected(u, v, c);

        for (int s = 1; s < n_; ++s) {
            const int t = parent_[static_cast<std::size_t>(s)];
            din.reset();
            const std::int64_t f = din.max_flow(s, t);
            const std::vector<char> side = din.min_cut_side(s);
            fl_[static_cast<std::size_t>(s)] = f;

            // every vertex still parented at t and sitting on s's side moves
            // under s. this is what splits the tree along the fresh cut.
            for (int v = 0; v < n_; ++v) {
                if (v != s && side[static_cast<std::size_t>(v)] &&
                    parent_[static_cast<std::size_t>(v)] == t) {
                    parent_[static_cast<std::size_t>(v)] = s;
                }
            }

            // grandparent rotation: t's parent on s's side means s is the higher
            // of the two, so lift s above t and swap the two edge weights.
            const int gp = parent_[static_cast<std::size_t>(t)];
            if (side[static_cast<std::size_t>(gp)]) {
                parent_[static_cast<std::size_t>(s)] = gp;
                parent_[static_cast<std::size_t>(t)] = s;
                fl_[static_cast<std::size_t>(s)] = fl_[static_cast<std::size_t>(t)];
                fl_[static_cast<std::size_t>(t)] = f;
            }
        }
    }

    // fill answer[u][v] = min tree-edge weight on the u..v path. n <= 150, so a
    // BFS from each source over the n-1 tree edges is O(n^2) total -- cheaper
    // than one max-flow, and it turns every query into an array read.
    void precompute_pairs() {
        // tree adjacency straight from the parent array: edge (i, parent[i])
        // with weight fl[i], for i = 1..n-1.
        std::vector<std::vector<std::pair<int, std::int64_t>>> adj(
            static_cast<std::size_t>(n_));
        for (int i = 1; i < n_; ++i) {
            const int p = parent_[static_cast<std::size_t>(i)];
            const std::int64_t w = fl_[static_cast<std::size_t>(i)];
            adj[static_cast<std::size_t>(i)].emplace_back(p, w);
            adj[static_cast<std::size_t>(p)].emplace_back(i, w);
        }

        answer_.assign(static_cast<std::size_t>(n_),
                       std::vector<std::int64_t>(static_cast<std::size_t>(n_), 0));
        std::vector<int> stack;
        for (int s = 0; s < n_; ++s) {
            std::vector<char> seen(static_cast<std::size_t>(n_), 0);
            std::vector<std::int64_t> best(static_cast<std::size_t>(n_), 0);
            seen[static_cast<std::size_t>(s)] = 1;
            stack.clear();
            stack.push_back(s);
            while (!stack.empty()) {
                const int v = stack.back();
                stack.pop_back();
                for (const auto& [u, w] : adj[static_cast<std::size_t>(v)]) {
                    if (seen[static_cast<std::size_t>(u)]) continue;
                    seen[static_cast<std::size_t>(u)] = 1;
                    // path-min: the bottleneck so far, tightened by this edge. v
                    // is either the source (no edge yet -> take w) or already has
                    // a running min.
                    const std::int64_t through =
                        (v == s) ? w
                                 : (best[static_cast<std::size_t>(v)] < w
                                        ? best[static_cast<std::size_t>(v)]
                                        : w);
                    best[static_cast<std::size_t>(u)] = through;
                    stack.push_back(u);
                }
            }
            answer_[static_cast<std::size_t>(s)] = best;
        }
    }

    int n_;
    std::vector<int> parent_;
    std::vector<std::int64_t> fl_;
    std::vector<std::vector<std::int64_t>> answer_;
};

}  // namespace p0038
