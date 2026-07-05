#pragma once

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <tuple>
#include <vector>

namespace p0032 {

// maximum weight matching in a GENERAL graph -- edmonds' primal-dual blossom,
// O(n^3). the matching need not be perfect; leaving a vertex unmatched is free,
// and the empty matching scores 0. weights are non-negative int; the answer and
// every dual live in int64 (n * 1e9 ~ 2.5e11 overflows int32, and a reduced cost
// lab[u]+lab[v]-2w reaches ~2e9 on its own).
//
// why this is the hard one -- bipartite matching has no odd cycles, so the
// hungarian method / min-cost-flow duals stay tight on a plain vertex cover.
// drop the bipartite promise and odd cycles appear. an augmenting search hits a
// vertex from two sides at once and forms an odd "blossom" -- a cycle you cannot
// two-color. flow can't model it. edmonds' answer: CONTRACT the blossom to a
// single pseudo-node, keep searching, and EXPAND it when the duals say so. that
// shrink/expand plus the dual bookkeeping is the whole difficulty.
//
// the primal-dual skeleton, all of it here:
//   dual -- every vertex carries a label lab[u]; every blossom a label lab[b].
//     feasibility is lab[u]+lab[v] >= 2*w(u,v) for every edge, tight edges
//     (reduced cost 0) are the only ones the matching may use. start every
//     vertex at w_max so all edges start slack.
//   grow -- from each exposed vertex plant an alternating tree over tight edges.
//     S=0 marks an outer node (even depth), S=1 an inner node (odd depth).
//   shrink -- a tight edge between two outer nodes with a common ancestor closes
//     an odd cycle. add_blossom contracts it; its odd length is what makes it a
//     blossom and not just a re-found path.
//   dual adjust -- no tight edge to grow on? slide the labels by the largest d
//     that keeps feasibility: outer labels down, inner up, blossom labels by 2d.
//     the moment an exposed outer real vertex would drop to <= 0 we STOP -- that
//     vertex stays unmatched. that single line is the non-perfect handling: a
//     zero label means "matching me further can't add weight," so we don't.
//   expand -- a blossom whose label hits 0 is torn back open, its internal
//     alternating structure spliced into the tree.
//   augment -- a tight edge between two outer nodes in DIFFERENT trees is an
//     augmenting path; flip it and the matching grows by one.
// each stage is O(n^2), at most O(n) stages -> O(n^3).
//
// sparsity -- we don't physically complete the graph. absent pairs sit at w=0
// and the scan skips every w<=0 edge, so a 0-weight (or missing) edge is never
// matched. that's identical to padding to a complete graph with 0-weight edges
// (a 0 edge never helps), just without the memory.
//
// LOWER BOUND -- reading the graph is Omega(n + m). the algorithm is O(n^3),
// which is the standard bound for dense weighted general matching; you can't go
// below reading the input, and no o(n^3) general-weighted-matching bound is
// known at this size. state it and move on.
//
// ML / systems note -- general-graph matching is the optimal pairing when the
// structure ISN'T bipartite: teammate pairing, sensor-to-sensor assignment,
// clustering into pairs. exactly the shape the clean flow reductions can't touch,
// because the odd cycle has no bipartite analog.
class Blossom {
public:
    using ll = std::int64_t;

    // build once from a 0-indexed edge list, dedup by max, then run.
    ll solve(int n, const std::vector<std::tuple<int, int, ll>>& edges) {
        n_ = n;
        if (n_ <= 0) return 0;
        // blossom pseudo-nodes live above n. at most one live pseudo-node per two
        // contracted base vertices, nested, so ids stay under 2n -- size for it.
        const int cap = 2 * n_ + 2;
        g_.assign(static_cast<std::size_t>(cap),
                  std::vector<Edge>(static_cast<std::size_t>(cap), Edge{0, 0, 0}));
        for (int u = 1; u <= n_; ++u)
            for (int v = 1; v <= n_; ++v) g_[u][v] = Edge{u, v, 0};

        ll w_max = 0;
        for (const auto& [a, b, w] : edges) {
            if (a == b || w <= 0) continue;  // no self-loops, ignore 0-weight.
            const int u = a + 1, v = b + 1;  // internal ids are 1-indexed.
            if (w > g_[u][v].w) {            // dedup parallel edges by max.
                g_[u][v].w = w;
                g_[v][u].w = w;
            }
            w_max = std::max(w_max, w);
        }

        lab_.assign(static_cast<std::size_t>(cap), 0);
        slack_.assign(static_cast<std::size_t>(cap), 0);
        match_.assign(static_cast<std::size_t>(cap), 0);
        st_.assign(static_cast<std::size_t>(cap), 0);
        pa_.assign(static_cast<std::size_t>(cap), 0);
        S_.assign(static_cast<std::size_t>(cap), 0);
        vis_.assign(static_cast<std::size_t>(cap), 0);
        flower_.assign(static_cast<std::size_t>(cap), {});
        flower_from_.assign(
            static_cast<std::size_t>(cap),
            std::vector<int>(static_cast<std::size_t>(cap), 0));

        n_x_ = n_;
        t_ = 0;
        for (int u = 0; u <= n_; ++u) {
            st_[u] = u;
            flower_[u].clear();
        }
        for (int u = 1; u <= n_; ++u)
            for (int v = 1; v <= n_; ++v)
                flower_from_[u][v] = (u == v ? u : 0);
        for (int u = 1; u <= n_; ++u) lab_[u] = w_max;

        while (run_stage()) {
            // each true is one augmentation -- another matched pair. false means
            // no profitable augmentation remains and the matching is optimal.
        }

        ll total = 0;
        for (int u = 1; u <= n_; ++u)
            if (match_[u] && match_[u] < u) total += g_[u][match_[u]].w;
        return total;
    }

private:
    struct Edge {
        int u, v;
        ll w;
    };

    static constexpr ll kInf = std::numeric_limits<ll>::max() / 4;

    int n_ = 0, n_x_ = 0, t_ = 0;
    std::vector<std::vector<Edge>> g_;
    std::vector<ll> lab_;     // dual labels on vertices and blossoms.
    std::vector<int> slack_;  // slack_[x] = the outer vertex giving x its tightest edge.
    std::vector<int> match_, st_, pa_, S_, vis_;
    std::vector<std::vector<int>> flower_from_;
    std::vector<std::vector<int>> flower_;
    std::deque<int> q_;

    // reduced cost of an edge under the current duals. feasibility keeps it >= 0;
    // a tight edge (== 0) is the only kind the matching is allowed to use.
    ll e_delta(const Edge& e) const {
        return lab_[e.u] + lab_[e.v] - g_[e.u][e.v].w * 2;
    }

    // slack_[x] tracks which outer real vertex u gives x its tightest edge.
    void update_slack(int u, int x) {
        if (!slack_[x] || e_delta(g_[u][x]) < e_delta(g_[slack_[x]][x]))
            slack_[x] = u;
    }
    void set_slack(int x) {
        slack_[x] = 0;
        for (int u = 1; u <= n_; ++u)
            if (g_[u][x].w > 0 && st_[u] != x && S_[st_[u]] == 0)
                update_slack(u, x);
    }

    // pushing a blossom onto the queue means pushing every base vertex under it.
    void q_push(int x) {
        if (x <= n_) {
            q_.push_back(x);
        } else {
            for (int c : flower_[x]) q_push(c);
        }
    }

    // relabel a subtree of base vertices to belong to blossom b.
    void set_st(int x, int b) {
        st_[x] = b;
        if (x > n_)
            for (int c : flower_[x]) set_st(c, b);
    }

    // position of xr inside b's cycle, rotated so the walk comes out even.
    int get_pr(int b, int xr) {
        int pr = static_cast<int>(
            std::find(flower_[b].begin(), flower_[b].end(), xr) -
            flower_[b].begin());
        if (pr % 2 == 1) {
            std::reverse(flower_[b].begin() + 1, flower_[b].end());
            return static_cast<int>(flower_[b].size()) - pr;
        }
        return pr;
    }

    // flip the match along the way into a blossom, then rotate it so its base is
    // the newly matched vertex.
    void set_match(int u, int v) {
        match_[u] = g_[u][v].v;
        if (u > n_) {
            Edge e = g_[u][v];
            int xr = flower_from_[u][e.u], pr = get_pr(u, xr);
            for (int i = 0; i < pr; ++i)
                set_match(flower_[u][static_cast<std::size_t>(i)],
                          flower_[u][static_cast<std::size_t>(i ^ 1)]);
            set_match(xr, v);
            std::rotate(flower_[u].begin(), flower_[u].begin() + pr,
                        flower_[u].end());
        }
    }

    // walk an augmenting path top-node to top-node, flipping matched edges.
    void augment(int u, int v) {
        for (;;) {
            int xnv = st_[match_[u]];
            set_match(u, v);
            if (!xnv) return;
            set_match(xnv, st_[pa_[xnv]]);
            u = st_[pa_[xnv]];
            v = xnv;
        }
    }

    // lowest common ancestor of two outer nodes in the alternating forest, or 0
    // if they sit in different trees (that's the augment case).
    int get_lca(int u, int v) {
        ++t_;
        for (; u || v; std::swap(u, v)) {
            if (u == 0) continue;
            if (vis_[u] == t_) return u;
            vis_[u] = t_;
            u = st_[match_[u]];
            if (u) u = st_[pa_[u]];
        }
        return 0;
    }

    // contract the odd cycle u..lca..v into one pseudo-node.
    void add_blossom(int u, int lca, int v) {
        int b = n_ + 1;
        while (b <= n_x_ && st_[b]) ++b;
        if (b > n_x_) ++n_x_;
        lab_[b] = 0;
        S_[b] = 0;
        match_[b] = match_[lca];
        flower_[b].clear();
        flower_[b].push_back(lca);
        for (int x = u, y; x != lca; x = st_[pa_[y]]) {
            flower_[b].push_back(x);
            y = st_[match_[x]];
            flower_[b].push_back(y);
            q_push(y);
        }
        std::reverse(flower_[b].begin() + 1, flower_[b].end());
        for (int x = v, y; x != lca; x = st_[pa_[y]]) {
            flower_[b].push_back(x);
            y = st_[match_[x]];
            flower_[b].push_back(y);
            q_push(y);
        }
        set_st(b, b);
        for (int x = 1; x <= n_x_; ++x) g_[b][x].w = g_[x][b].w = 0;
        for (int x = 1; x <= n_; ++x) flower_from_[b][x] = 0;
        for (int xs : flower_[b]) {
            for (int x = 1; x <= n_x_; ++x)
                if (g_[b][x].w == 0 || e_delta(g_[xs][x]) < e_delta(g_[b][x])) {
                    g_[b][x] = g_[xs][x];
                    g_[x][b] = g_[x][xs];
                }
            for (int x = 1; x <= n_; ++x)
                if (flower_from_[xs][x]) flower_from_[b][x] = xs;
        }
        set_slack(b);
    }

    // tear a blossom back open, splicing its alternating chain into the tree.
    void expand_blossom(int b) {
        for (int c : flower_[b]) set_st(c, c);
        int xr = flower_from_[b][g_[b][pa_[b]].u], pr = get_pr(b, xr);
        for (int i = 0; i < pr; i += 2) {
            int xs = flower_[b][static_cast<std::size_t>(i)];
            int xns = flower_[b][static_cast<std::size_t>(i + 1)];
            pa_[xs] = g_[xns][xs].u;
            S_[xs] = 1;
            S_[xns] = 0;
            slack_[xs] = 0;
            set_slack(xns);
            q_push(xns);
        }
        S_[xr] = 1;
        pa_[xr] = pa_[b];
        for (std::size_t i = static_cast<std::size_t>(pr) + 1;
             i < flower_[b].size(); ++i) {
            int xs = flower_[b][i];
            S_[xs] = -1;
            set_slack(xs);
        }
        st_[b] = 0;
    }

    // a tight edge found during a grow. extend the tree, augment, or shrink.
    bool on_found_edge(const Edge& e) {
        int u = st_[e.u], v = st_[e.v];
        if (S_[v] == -1) {
            pa_[v] = e.u;
            S_[v] = 1;
            int nu = st_[match_[v]];
            slack_[v] = slack_[nu] = 0;
            S_[nu] = 0;
            q_push(nu);
        } else if (S_[v] == 0) {
            int lca = get_lca(u, v);
            if (!lca) {
                augment(u, v);
                augment(v, u);
                return true;  // matching grew -- restart the search.
            }
            add_blossom(u, lca, v);
        }
        return false;
    }

    // one primal-dual stage: grow trees from every exposed vertex, adjust duals
    // when stuck, and return true the instant an augmenting path appears. false
    // means the current matching is already maximum-weight.
    bool run_stage() {
        for (int i = 1; i <= n_x_; ++i) {
            S_[i] = -1;
            slack_[i] = 0;
        }
        q_.clear();
        for (int x = 1; x <= n_x_; ++x)
            if (st_[x] == x && !match_[x]) {
                pa_[x] = 0;
                S_[x] = 0;
                q_push(x);
            }
        if (q_.empty()) return false;

        for (;;) {
            while (!q_.empty()) {
                int u = q_.front();
                q_.pop_front();
                if (S_[st_[u]] == 1) continue;
                for (int v = 1; v <= n_; ++v)
                    if (g_[u][v].w > 0 && st_[u] != st_[v]) {
                        if (e_delta(g_[u][v]) == 0) {
                            if (on_found_edge(g_[u][v])) return true;
                        } else {
                            update_slack(u, st_[v]);
                        }
                    }
            }

            // no tight edge to grow on -- slide the duals by the largest safe d.
            ll d = kInf;
            for (int b = n_ + 1; b <= n_x_; ++b)
                if (st_[b] == b && S_[b] == 1) d = std::min(d, lab_[b] / 2);
            for (int x = 1; x <= n_x_; ++x)
                if (st_[x] == x && slack_[x]) {
                    if (S_[x] == -1)
                        d = std::min(d, e_delta(g_[slack_[x]][x]));
                    else if (S_[x] == 0)
                        d = std::min(d, e_delta(g_[slack_[x]][x]) / 2);
                }
            for (int u = 1; u <= n_; ++u) {
                if (S_[st_[u]] == 0) {
                    if (lab_[u] <= d) return false;  // <-- non-perfect stop.
                    lab_[u] -= d;
                } else if (S_[st_[u]] == 1) {
                    lab_[u] += d;
                }
            }
            for (int b = n_ + 1; b <= n_x_; ++b)
                if (st_[b] == b) {
                    if (S_[b] == 0)
                        lab_[b] += d * 2;
                    else if (S_[b] == 1)
                        lab_[b] -= d * 2;
                }

            q_.clear();
            for (int x = 1; x <= n_x_; ++x)
                if (st_[x] == x && slack_[x] && st_[slack_[x]] != x &&
                    e_delta(g_[slack_[x]][x]) == 0)
                    if (on_found_edge(g_[slack_[x]][x])) return true;
            for (int b = n_ + 1; b <= n_x_; ++b)
                if (st_[b] == b && S_[b] == 1 && lab_[b] == 0)
                    expand_blossom(b);
        }
    }
};

// one-shot entry. 0-indexed edges (u, v, w); parallel edges deduped by max,
// self-loops and non-positive weights ignored. returns the max matching weight.
inline std::int64_t max_weight_matching(
    int n, const std::vector<std::tuple<int, int, std::int64_t>>& edges) {
    Blossom b;
    return b.solve(n, edges);
}

}  // namespace p0032
