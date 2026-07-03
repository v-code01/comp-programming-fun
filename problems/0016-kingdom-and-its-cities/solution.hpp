#pragma once

#include <algorithm>
#include <utility>
#include <vector>

namespace p0016 {

// min cities to capture so no two important cities stay connected -- per query,
// the important set changes but the tree does not. the naive move is "run a cut
// dp over all n nodes for every query": that's O(n) per query and sum k up to
// 1e5 queries blows it to 1e10. the whole trick is that a query only cares about
// the important nodes and the junctions where their paths meet -- everything else
// on a path is an interchangeable non-important node you can cut. compress the
// tree down to those O(k) relevant nodes -- the VIRTUAL (auxiliary) tree -- and
// dp on that instead.
//
// the -1 rule: you may only capture NON-important cities, so two important cities
// joined by a direct edge can never be separated -- nothing sits between them to
// remove. rooted, an edge is a parent link, so that's exactly "some important
// node's parent is also important". catch it in O(k) and bail.
//
// the dp, bottom-up on the virtual tree, returns per node whether its subtree
// still has an "exposed" important node reaching up through it:
//   v important        -> every exposed child path connects to v and must be cut,
//                         so answer += (#exposed children); v itself stays exposed.
//   v not important,
//     >=2 exposed kids -> cut v once, sever all of them at the junction; not exposed.
//     ==1 exposed kid  -> pass it upward, defer -- a higher cut may serve several
//                         branches with one node; still exposed.
//     ==0              -> nothing to do; not exposed.
// deferring is why the single-child case passes up instead of cutting: the one
// exposed chain either meets an important ancestor (cut once there anyway) or a
// >=2 junction (one cut severs many) -- never worse than cutting eagerly.
//
// build tree + lca tables once: O(n log n). per query: sort k by dfs-time
// O(k log k), k-1 lca calls O(k log n), stack-build + dp + cleanup O(k). total
// O(n log n + sum k * log n). the log is the lca jump; drop a better lca (euler +
// O(1) sparse table) and per query is O(k log k) for the sort alone.
//
// lower bound: reading the tree and every query is Omega(n + sum k). per query
// you must at least sort the k important nodes to know their meeting structure,
// Omega(k log k) comparison-model, so sum is Omega(n + sum k log k). this hits it
// up to the lca log factor.
class Kingdom {
public:
    // adjacency + rooted structure + binary-lifting lca. nodes are 1..n. edges
    // are 0-indexed pairs of node ids. O(n log n) time and space.
    void build(int n, const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        adj_.assign(static_cast<std::size_t>(n) + 1, {});
        for (const auto& e : edges) {
            adj_[static_cast<std::size_t>(e.first)].push_back(e.second);
            adj_[static_cast<std::size_t>(e.second)].push_back(e.first);
        }

        // one past the top bit of n -- enough ancestors to jump any depth.
        log_ = 1;
        while ((1 << log_) < n_) ++log_;
        ++log_;

        tin_.assign(static_cast<std::size_t>(n) + 1, 0);
        tout_.assign(static_cast<std::size_t>(n) + 1, 0);
        depth_.assign(static_cast<std::size_t>(n) + 1, 0);
        up_.assign(static_cast<std::size_t>(log_),
                   std::vector<int>(static_cast<std::size_t>(n) + 1, kRoot));

        dfs_iterative();

        // up_[j][v] = the 2^j-th ancestor. root points at itself so jumps past
        // the top saturate instead of walking off the array.
        for (int j = 1; j < log_; ++j) {
            const auto& prev = up_[static_cast<std::size_t>(j - 1)];
            auto& cur = up_[static_cast<std::size_t>(j)];
            for (int v = 1; v <= n_; ++v) {
                cur[static_cast<std::size_t>(v)] =
                    prev[static_cast<std::size_t>(prev[static_cast<std::size_t>(v)])];
            }
        }

        // per-query scratch, all n-sized but only ever touched on O(k) entries and
        // reset on those same entries -- never an O(n) wipe inside a query.
        is_imp_.assign(static_cast<std::size_t>(n) + 1, 0);
        cnt_.assign(static_cast<std::size_t>(n) + 1, 0);
        vpar_.assign(static_cast<std::size_t>(n) + 1, 0);
    }

    // answer one query. imp is the important set; it gets sorted in place. returns
    // the min captures, or -1 if two important cities are adjacent.
    int query(std::vector<int>& imp) {
        const int k = static_cast<int>(imp.size());
        if (k <= 1) return 0;  // one (or no) important city -- nothing to separate.

        for (int a : imp) is_imp_[static_cast<std::size_t>(a)] = 1;

        // adjacency -> impossible. an edge is a parent link once rooted, so check
        // each important node's parent. do it before the build so we can bail and
        // still clean up only the k marks.
        int impossible = 0;
        for (int a : imp) {
            if (a != kRoot &&
                is_imp_[static_cast<std::size_t>(up_[0][static_cast<std::size_t>(a)])]) {
                impossible = 1;
                break;
            }
        }
        if (impossible) {
            for (int a : imp) is_imp_[static_cast<std::size_t>(a)] = 0;
            return -1;
        }

        // ---- assemble the virtual tree's node set: important + pairwise-meet
        // lcas. sorting by dfs entry time puts a node right before its whole
        // subtree, so consecutive pairs are enough -- their lcas close the set
        // under lca (standard auxiliary-tree fact). ----
        std::sort(imp.begin(), imp.end(), [this](int x, int y) {
            return tin_[static_cast<std::size_t>(x)] < tin_[static_cast<std::size_t>(y)];
        });

        vnodes_.clear();
        for (int a : imp) vnodes_.push_back(a);
        for (int i = 0; i + 1 < k; ++i) vnodes_.push_back(lca(imp[i], imp[i + 1]));

        std::sort(vnodes_.begin(), vnodes_.end(), [this](int x, int y) {
            return tin_[static_cast<std::size_t>(x)] < tin_[static_cast<std::size_t>(y)];
        });
        // unique by node id -- tin is a bijection, so equal tin means same node.
        vnodes_.erase(std::unique(vnodes_.begin(), vnodes_.end(),
                                  [this](int x, int y) {
                                      return tin_[static_cast<std::size_t>(x)] ==
                                             tin_[static_cast<std::size_t>(y)];
                                  }),
                      vnodes_.end());

        // ---- stack-build parents. walk nodes in tin order; the stack holds the
        // current ancestor chain. pop until the top is an ancestor of u -- that
        // top is u's virtual parent. nodes[0] is the min-tin node, ancestor of
        // all, so the stack never empties under it. ----
        const int m = static_cast<int>(vnodes_.size());
        stack_.clear();
        stack_.push_back(vnodes_[0]);
        vpar_[static_cast<std::size_t>(vnodes_[0])] = kNone;  // the virtual root.
        for (int i = 1; i < m; ++i) {
            const int u = vnodes_[static_cast<std::size_t>(i)];
            while (!is_ancestor(stack_.back(), u)) stack_.pop_back();
            vpar_[static_cast<std::size_t>(u)] = stack_.back();
            stack_.push_back(u);
        }

        // ---- the dp. reverse tin order visits every child before its parent
        // (a parent's tin is strictly smaller), so cnt_[v] is fully accumulated by
        // the time we read it. ----
        int answer = 0;
        for (int i = m - 1; i >= 0; --i) {
            const int v = vnodes_[static_cast<std::size_t>(i)];
            const int exposed = cnt_[static_cast<std::size_t>(v)];
            int up_here;
            if (is_imp_[static_cast<std::size_t>(v)]) {
                answer += exposed;  // one cut per exposed child path into v.
                up_here = 1;        // v is important and uncut -- reaches upward.
            } else if (exposed >= 2) {
                ++answer;           // cut v, sever every child branch at once.
                up_here = 0;
            } else {
                up_here = exposed;  // 1 -> pass up (defer); 0 -> nothing.
            }
            const int p = vpar_[static_cast<std::size_t>(v)];
            if (up_here && p != kNone) ++cnt_[static_cast<std::size_t>(p)];
        }

        // ---- cleanup: reset only the entries we touched. is_imp_ on the k
        // important, cnt_ on the m virtual nodes. O(k), never O(n). ----
        for (int a : imp) is_imp_[static_cast<std::size_t>(a)] = 0;
        for (int v : vnodes_) cnt_[static_cast<std::size_t>(v)] = 0;

        return answer;
    }

private:
    static constexpr int kRoot = 1;   // the tree is rooted here.
    static constexpr int kNone = 0;   // "no parent" sentinel (node 0 is unused).

    // iterative dfs -- n up to 1e5 would blow a recursive stack. assigns entry and
    // exit stamps from one running timer plus depth and the immediate parent.
    void dfs_iterative() {
        std::vector<int> it(static_cast<std::size_t>(n_) + 1, 0);  // adj cursor per node.
        std::vector<int> st;
        st.reserve(static_cast<std::size_t>(n_));
        int timer = 0;

        up_[0][kRoot] = kRoot;  // root's parent is itself.
        depth_[static_cast<std::size_t>(kRoot)] = 0;
        tin_[static_cast<std::size_t>(kRoot)] = timer++;
        st.push_back(kRoot);

        while (!st.empty()) {
            const int v = st.back();
            auto& cursor = it[static_cast<std::size_t>(v)];
            const auto& nbrs = adj_[static_cast<std::size_t>(v)];
            if (cursor < static_cast<int>(nbrs.size())) {
                const int u = nbrs[static_cast<std::size_t>(cursor++)];
                if (u == up_[0][static_cast<std::size_t>(v)] && v != kRoot) continue;
                if (u == v) continue;  // paranoia: no self-loops in a tree.
                up_[0][static_cast<std::size_t>(u)] = v;
                depth_[static_cast<std::size_t>(u)] =
                    depth_[static_cast<std::size_t>(v)] + 1;
                tin_[static_cast<std::size_t>(u)] = timer++;
                st.push_back(u);
            } else {
                tout_[static_cast<std::size_t>(v)] = timer++;
                st.pop_back();
            }
        }
    }

    // a is an ancestor of b iff b's whole [tin,tout] window nests inside a's.
    bool is_ancestor(int a, int b) const {
        return tin_[static_cast<std::size_t>(a)] <= tin_[static_cast<std::size_t>(b)] &&
               tout_[static_cast<std::size_t>(b)] <= tout_[static_cast<std::size_t>(a)];
    }

    int lca(int u, int v) const {
        if (depth_[static_cast<std::size_t>(u)] < depth_[static_cast<std::size_t>(v)]) {
            std::swap(u, v);
        }
        int diff = depth_[static_cast<std::size_t>(u)] - depth_[static_cast<std::size_t>(v)];
        for (int j = 0; j < log_; ++j) {
            if ((diff >> j) & 1) u = up_[static_cast<std::size_t>(j)][static_cast<std::size_t>(u)];
        }
        if (u == v) return u;
        for (int j = log_ - 1; j >= 0; --j) {
            const int au = up_[static_cast<std::size_t>(j)][static_cast<std::size_t>(u)];
            const int av = up_[static_cast<std::size_t>(j)][static_cast<std::size_t>(v)];
            if (au != av) {
                u = au;
                v = av;
            }
        }
        return up_[0][static_cast<std::size_t>(u)];
    }

    int n_ = 0;
    int log_ = 1;
    std::vector<std::vector<int>> adj_;
    std::vector<int> tin_, tout_, depth_;
    std::vector<std::vector<int>> up_;  // up_[j][v] = 2^j-th ancestor of v.

    // per-query scratch, reset on touched entries only.
    std::vector<int> is_imp_;
    std::vector<int> cnt_;
    std::vector<int> vpar_;
    std::vector<int> vnodes_;  // the virtual tree's nodes, tin-sorted.
    std::vector<int> stack_;   // ancestor chain during the stack build.
};

// convenience for tests: build once, answer a batch of queries. each query vector
// is taken by value so callers keep their originals.
inline std::vector<int> solve(int n, const std::vector<std::pair<int, int>>& edges,
                              std::vector<std::vector<int>> queries) {
    Kingdom k;
    k.build(n, edges);
    std::vector<int> out;
    out.reserve(queries.size());
    for (auto& q : queries) out.push_back(k.query(q));
    return out;
}

}  // namespace p0016
