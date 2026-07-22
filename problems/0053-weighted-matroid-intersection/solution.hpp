#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace p0053 {

// maximum-weight colorful spanning forest -- WEIGHTED MATROID INTERSECTION.
//
// same two matroids as the unweighted twin (problem 0048), same m edges:
//   (A) forest -- no cycle. the GRAPHIC matroid M1.
//   (B) colorful -- at most one edge per color. a PARTITION matroid M2.
// but now every edge carries a weight, any integer, negatives allowed. pick the
// edge set independent in BOTH that MAXIMIZES total weight. the set need not be
// large -- if every weight is negative the empty set (weight 0) wins outright.
//
// the unweighted version just wanted the biggest common independent set, and a
// SHORTEST augmenting path delivered it. weight breaks that. you no longer want
// the most edges -- you want the heaviest legal pile, and a heavy edge three
// swaps away can beat a light one you could add for free. so the path rule
// changes: not fewest arcs, but MAX PROFIT, and fewest arcs only to break ties.
// that swap is the whole difficulty of the weighted half.
//
// W(k) -- the max weight over all common independent sets of size EXACTLY k --
// is CONCAVE in k. that single fact runs the whole method:
//   W(0) = 0. the empty set is optimal at size 0, trivially.
//   given the optimal size-k set I, the best size-(k+1) set is I augmented along
//   the MAX-PROFIT source->sink path of the exchange graph. one path, one edge
//   more, and I stays optimal for its new size.
//   concavity means the marginal gain W(k+1)-W(k) only shrinks. so the first
//   time the best path's profit drops to <= 0, every later one is <= it too --
//   nothing ahead can lift the total. stop. the peak is the answer.
//
// the exchange graph D on the m edges, rebuilt around I every step (the GATING
// is identical to 0048 -- what changes is how we walk it):
//   sources  X1 = { y not in I : I+y is a forest }        -- free to add in M1.
//   sinks    X2 = { y not in I : color(y) unused in I }    -- free to add in M2.
//   arc x->y  (x in I, y not in I)  iff  I-x+y is a forest    -- M1 gates it.
//   arc y->x  (y not in I, x in I)  iff  I-x+y is colorful    -- M2 gates it.
// give each element a PROFIT: +w(y) for an out-node y (you'd ADD it), -w(x) for
// an in-node x (you'd DROP it). a source->sink path's profit is the sum over its
// nodes, and it equals W(new) - W(old) exactly -- add the out-nodes, drop the
// in-nodes, |I| climbs by one. maximize that.
//
// MAX-PROFIT, then FEWEST ARCS. the tie-break is load-bearing, not cosmetic. at
// the optimum the exchange graph has NO positive-profit cycle (that's the
// invariant "I is max-weight for its size" made geometric). a zero-profit cycle
// can still exist -- and routing through one would corrupt I into a dependent
// set without changing the weight. fewest-arcs forbids it: a cycle only adds
// arcs, never removes them, so the minimum-arc path never loops. drop this
// tie-break and the solver silently returns dependent sets. this is the bug.
//
// no positive cycle + a strict penalty on extra arcs => the "longest path" here
// is well posed. BELLMAN-FORD relaxes it: dist keyed by (profit desc, arcs asc),
// seeded at every source with its own profit, read off at the best sink. no
// non-negative cycle, so it converges in at most one pass per path edge.
//
// NOTE -- the unweighted twin (0048) could PRUNE the exchange graph: skip arcs
// into sources, skip a sink's out-arcs, because a SHORTEST path never uses them.
// here it can. a max-profit path may thread through a source or past a sink to
// reach heavier weight, so we build the FULL arc set. that density is the price
// of weights.
//
// COMPLEXITY -- one augmentation: an O(m^2) arc build (a source contributes |I|
// arcs; a fundamental cycle scan is O(n) per non-source) plus Bellman-Ford over
// E = O(m^2) arcs for O(m) passes -> O(m^3) worst. augmentations <= r =
// min(n-1, #colors) <= m, and concavity stops us early. total O(r * m^3), and at
// n=m=200 the early stop and short paths keep it far under the worst case.
//
// LOWER BOUND -- reading the graph is Omega(n+m); no answer without touching
// every edge. above that, weighted matroid intersection is polynomial. the
// concavity of W(k) is the structural heart: it both proves the augmenting
// scheme optimal at each size AND lets you halt the instant W stops rising.
// unweighted intersection needs only shortest augmenting paths; the weights are
// what force max-profit-then-fewest-arcs, the harder half. three matroids would
// be NP-hard -- two is the tractable frontier.
class WeightedMatroidIntersection {
public:
    using ll = std::int64_t;
    // (u, v, color, weight); u,v in 1..n, color >= 1, weight any int64.
    using Edge = std::tuple<int, int, int, ll>;

    // the max achievable total weight. always >= 0 -- the empty set is legal.
    ll solve(int n, const std::vector<Edge>& edges) {
        n_ = n;
        m_ = static_cast<int>(edges.size());
        if (m_ == 0 || n_ <= 0) return 0;

        eu_.resize(static_cast<std::size_t>(m_));
        ev_.resize(static_cast<std::size_t>(m_));
        ec_.resize(static_cast<std::size_t>(m_));
        ew_.resize(static_cast<std::size_t>(m_));
        int max_color = 1;
        for (int j = 0; j < m_; ++j) {
            const auto& [u, v, c, w] = edges[static_cast<std::size_t>(j)];
            eu_[static_cast<std::size_t>(j)] = u - 1;  // internal 0-indexed.
            ev_[static_cast<std::size_t>(j)] = v - 1;
            ec_[static_cast<std::size_t>(j)] = c;
            ew_[static_cast<std::size_t>(j)] = w;
            if (c > max_color) max_color = c;
        }
        max_color_ = max_color;

        in_i_.assign(static_cast<std::size_t>(m_), 0);
        ll cur = 0;    // W(|I|) -- weight of the current optimal-for-its-size set.
        ll best = 0;   // the empty set is always on the table.
        std::vector<int> path;
        ll profit = 0;
        while (find_augment(path, profit)) {
            // W is concave -- once the best marginal gain is non-positive, every
            // later one is too. nothing ahead beats what we already have. stop.
            if (profit <= 0) break;
            for (const int v : path) in_i_[static_cast<std::size_t>(v)] ^= 1;
            cur += profit;
            if (cur > best) best = cur;
        }
        return best;
    }

private:
    int n_ = 0, m_ = 0, max_color_ = 1;
    std::vector<int> eu_, ev_, ec_;  // endpoints (0-indexed) and colors.
    std::vector<ll> ew_;             // weights.
    std::vector<int> in_i_;          // membership of each edge in I.
    std::vector<int> dsu_;           // rebuilt every step over the edges in I.

    // path-halving find over the forest I. tiny structure, correctness first.
    int find(int x) {
        while (dsu_[static_cast<std::size_t>(x)] != x) {
            dsu_[static_cast<std::size_t>(x)] =
                dsu_[static_cast<std::size_t>(dsu_[static_cast<std::size_t>(x)])];
            x = dsu_[static_cast<std::size_t>(x)];
        }
        return x;
    }

    // build the exchange graph around I, find the MAX-PROFIT (tie: fewest-arcs)
    // source->sink path, hand it back. profit_out carries its weight delta. true
    // iff any source->sink path exists (I is not yet maximum cardinality).
    bool find_augment(std::vector<int>& path, ll& profit_out) {
        // --- DSU + forest adjacency + color holders over the current I. ---
        dsu_.assign(static_cast<std::size_t>(n_), 0);
        for (int i = 0; i < n_; ++i) dsu_[static_cast<std::size_t>(i)] = i;
        std::vector<std::vector<std::pair<int, int>>> adj(
            static_cast<std::size_t>(n_));  // node -> (neighbor, edge index).
        std::vector<int> color_holder(static_cast<std::size_t>(max_color_ + 1), -1);
        std::vector<int> i_set;  // the edges currently in I.
        for (int j = 0; j < m_; ++j) {
            if (!in_i_[static_cast<std::size_t>(j)]) continue;
            const int u = eu_[static_cast<std::size_t>(j)];
            const int v = ev_[static_cast<std::size_t>(j)];
            const int a = find(u), b = find(v);
            if (a != b) dsu_[static_cast<std::size_t>(a)] = b;  // I is a forest.
            adj[static_cast<std::size_t>(u)].push_back({v, j});
            adj[static_cast<std::size_t>(v)].push_back({u, j});
            color_holder[static_cast<std::size_t>(ec_[static_cast<std::size_t>(j)])] = j;
            i_set.push_back(j);
        }

        // --- the FULL arc set, keyed by tail. no shortest-path pruning here. ---
        std::vector<std::vector<int>> out(static_cast<std::size_t>(m_));
        // reusable stamped BFS buffers for the forest tree-path scans.
        std::vector<int> vis(static_cast<std::size_t>(n_), -1);
        std::vector<int> pv(static_cast<std::size_t>(n_), -1);
        std::vector<int> pe(static_cast<std::size_t>(n_), -1);
        int stamp = 0;
        std::vector<int> bq;

        for (int y = 0; y < m_; ++y) {
            if (in_i_[static_cast<std::size_t>(y)]) continue;  // out-nodes only.
            const int uy = eu_[static_cast<std::size_t>(y)];
            const int vy = ev_[static_cast<std::size_t>(y)];
            const bool is_source = (find(uy) != find(vy));
            const int cy = ec_[static_cast<std::size_t>(y)];
            const bool is_sink = (color_holder[static_cast<std::size_t>(cy)] == -1);

            // M2 arcs y->x: which x make I-x+y colorful again.
            if (is_sink) {
                // color(y) unused -> I+y already colorful -> I-x+y colorful for
                // EVERY x. a max-profit path may pass through y, so keep them all.
                for (const int x : i_set) out[static_cast<std::size_t>(y)].push_back(x);
            } else {
                // color(y) used by exactly one holder -- only dropping it helps.
                out[static_cast<std::size_t>(y)].push_back(
                    color_holder[static_cast<std::size_t>(cy)]);
            }

            // M1 arcs x->y: which x make I-x+y a forest.
            if (is_source) {
                // I+y is a forest -> I-x+y subset of it, a forest for EVERY x.
                for (const int x : i_set) out[static_cast<std::size_t>(x)].push_back(y);
            } else {
                // I+y closes one fundamental cycle = the tree path between y's
                // endpoints. removing an x on that path reopens it -> arc x->y.
                ++stamp;
                vis[static_cast<std::size_t>(uy)] = stamp;
                pv[static_cast<std::size_t>(uy)] = -1;
                pe[static_cast<std::size_t>(uy)] = -1;
                bq.clear();
                bq.push_back(uy);
                for (std::size_t h = 0; h < bq.size(); ++h) {
                    const int x = bq[h];
                    if (x == vy) break;
                    for (const auto& [nb, ei] : adj[static_cast<std::size_t>(x)]) {
                        if (vis[static_cast<std::size_t>(nb)] == stamp) continue;
                        vis[static_cast<std::size_t>(nb)] = stamp;
                        pv[static_cast<std::size_t>(nb)] = x;
                        pe[static_cast<std::size_t>(nb)] = ei;
                        bq.push_back(nb);
                    }
                }
                for (int cur = vy; cur != uy;
                     cur = pv[static_cast<std::size_t>(cur)]) {
                    const int x = pe[static_cast<std::size_t>(cur)];
                    out[static_cast<std::size_t>(x)].push_back(y);
                }
            }
        }

        // --- node profits: +w to add an out-node, -w to drop an in-node. ---
        std::vector<ll> prof(static_cast<std::size_t>(m_));
        for (int j = 0; j < m_; ++j)
            prof[static_cast<std::size_t>(j)] = in_i_[static_cast<std::size_t>(j)]
                                                    ? -ew_[static_cast<std::size_t>(j)]
                                                    : ew_[static_cast<std::size_t>(j)];

        // --- Bellman-Ford: maximize (profit, then minimize arcs). ---
        constexpr ll kNeg = std::numeric_limits<ll>::min() / 4;
        std::vector<ll> dprofit(static_cast<std::size_t>(m_), kNeg);
        std::vector<int> darcs(static_cast<std::size_t>(m_), 0);
        std::vector<int> par(static_cast<std::size_t>(m_), -1);

        // seed every source with its own profit at zero arcs. relaxation may
        // still improve a source if a heavier path happens to reach it.
        for (int y = 0; y < m_; ++y) {
            if (in_i_[static_cast<std::size_t>(y)]) continue;
            if (find(eu_[static_cast<std::size_t>(y)]) !=
                find(ev_[static_cast<std::size_t>(y)])) {
                dprofit[static_cast<std::size_t>(y)] = prof[static_cast<std::size_t>(y)];
                darcs[static_cast<std::size_t>(y)] = 0;
                par[static_cast<std::size_t>(y)] = -1;
            }
        }

        // no non-negative cycle -> a path never gains from a loop -> converges in
        // <= m passes. break the instant a full pass changes nothing.
        for (int pass = 0; pass <= m_; ++pass) {
            bool changed = false;
            for (int t = 0; t < m_; ++t) {
                if (dprofit[static_cast<std::size_t>(t)] == kNeg) continue;
                for (const int h : out[static_cast<std::size_t>(t)]) {
                    const ll np = dprofit[static_cast<std::size_t>(t)] +
                                  prof[static_cast<std::size_t>(h)];
                    const int na = darcs[static_cast<std::size_t>(t)] + 1;
                    ll& hp = dprofit[static_cast<std::size_t>(h)];
                    int& ha = darcs[static_cast<std::size_t>(h)];
                    if (np > hp || (np == hp && na < ha)) {
                        hp = np;
                        ha = na;
                        par[static_cast<std::size_t>(h)] = t;
                        changed = true;
                    }
                }
            }
            if (!changed) break;
        }

        // --- read off the best sink: max profit, fewest arcs to break ties. ---
        ll best_p = kNeg;
        int best_arcs = 0, target = -1;
        for (int y = 0; y < m_; ++y) {
            if (in_i_[static_cast<std::size_t>(y)]) continue;
            if (color_holder[static_cast<std::size_t>(ec_[static_cast<std::size_t>(y)])] !=
                -1)
                continue;  // not a sink.
            if (dprofit[static_cast<std::size_t>(y)] == kNeg) continue;
            const ll p = dprofit[static_cast<std::size_t>(y)];
            const int a = darcs[static_cast<std::size_t>(y)];
            if (p > best_p || (p == best_p && a < best_arcs)) {
                best_p = p;
                best_arcs = a;
                target = y;
            }
        }

        if (target < 0) return false;  // no source->sink path -- I is maximal.
        // walk parents back to the source; XOR flips membership along the path.
        path.clear();
        for (int cur = target; cur != -1; cur = par[static_cast<std::size_t>(cur)])
            path.push_back(cur);
        profit_out = best_p;
        return true;
    }
};

// one-shot entry. edges are (u, v, color, weight) with u,v in 1..n, color >= 1,
// weight any int64. parallel edges, multi-edges and self-loops are all accepted
// (a self-loop never joins a forest, so it never enters I). returns the maximum
// total weight of an edge set that is a forest AND uses each color at most once,
// which is >= 0 since the empty set is always allowed.
inline std::int64_t max_weight_colorful_forest(
    int n, const std::vector<std::tuple<int, int, int, std::int64_t>>& edges) {
    WeightedMatroidIntersection wmi;
    return wmi.solve(n, edges);
}

}  // namespace p0053
