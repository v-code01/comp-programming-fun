#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "solution.hpp"  // for p0033::Edge and p0033::i64 -- the POD types only.

namespace p0033 {

// the oracle. no lagrangian, no toll, no binary search -- just the definition of
// g(t), walked bottom-up as a subtree knapsack. for EVERY k at once.
//
// f[v] is two polynomials over "edges used inside subtree(v)":
//     f0[v][k] = max weight of a k-edge matching in subtree(v), v NOT matched.
//     f1[v][k] = max weight of a k-edge matching in subtree(v), v matched to a
//                child. (v may still get matched to its parent higher up -- that
//                decision is made where the parent folds v in.)
// fold each child c into v with the standard size-bounded knapsack:
//   * skip edge (v,c): add max(f0[c][b], f1[c][b]) at column a+b.
//   * take edge (v,c): needs v free (f0 side) and c free (f0[c]); lands at column
//     a+b+1 with +w. v becomes matched.
// merging two polynomials of widths |Tv|, |Tc| is O(|Tv|*|Tc|); summed over the
// tree that telescopes to O(n^2). sane for small n -- which is exactly where the
// differential test against the aliens solution lives.
//
// iterative, deepest-child-first, no recursion -- so even a 1e5 bamboo can't blow
// the stack, though the O(n^2) fold keeps this to small n in practice.
//
// returns ans[k] for k in [0, n/2]: the exact g(k), or -1 where no k-edge
// matching exists. weights are >= 1, so a feasible k >= 1 is always > 0 and -1
// can't collide with a real answer; k=0 is 0.
class TreeMatchRef {
public:
    TreeMatchRef(int n, const std::vector<Edge>& edges) : n_(n < 0 ? 0 : n) {
        build(edges);
        compute();
    }

    // exact g(k), or -1 if a k-edge matching can't exist.
    i64 at(int k) const {
        if (k < 0) return -1;
        if (k == 0) return 0;
        if (n_ <= 0) return -1;
        const std::vector<i64>& r0 = f0_[static_cast<std::size_t>(root_)];
        const std::vector<i64>& r1 = f1_[static_cast<std::size_t>(root_)];
        i64 best = NEG;
        if (static_cast<std::size_t>(k) < r0.size()) best = std::max(best, r0[static_cast<std::size_t>(k)]);
        if (static_cast<std::size_t>(k) < r1.size()) best = std::max(best, r1[static_cast<std::size_t>(k)]);
        return best <= NEG ? -1 : best;
    }

    // g(k) for every k in [0, n/2].
    std::vector<i64> all() const {
        std::vector<i64> ans(static_cast<std::size_t>(n_ / 2 + 1));
        for (std::size_t k = 0; k < ans.size(); ++k) ans[k] = at(static_cast<int>(k));
        return ans;
    }

private:
    static constexpr i64 NEG = -(static_cast<i64>(1) << 62);  // unreachable column.

    void build(const std::vector<Edge>& edges) {
        const std::size_t N = static_cast<std::size_t>(n_);
        parent_.assign(N, -1);
        wpar_.assign(N, 0);
        order_.clear();
        root_ = 0;
        if (n_ == 0) return;

        std::vector<int> deg(N, 0);
        for (const Edge& e : edges) {
            ++deg[static_cast<std::size_t>(e.u - 1)];
            ++deg[static_cast<std::size_t>(e.v - 1)];
        }
        std::vector<int> head(N + 1, 0);
        for (std::size_t v = 0; v < N; ++v) head[v + 1] = head[v] + deg[v];
        std::vector<int> adj(static_cast<std::size_t>(head[N]));
        std::vector<i64> aw(static_cast<std::size_t>(head[N]));
        std::vector<int> cur = head;
        for (const Edge& e : edges) {
            const int a = e.u - 1, b = e.v - 1;
            const std::size_t ia = static_cast<std::size_t>(cur[static_cast<std::size_t>(a)]++);
            const std::size_t ib = static_cast<std::size_t>(cur[static_cast<std::size_t>(b)]++);
            adj[ia] = b; aw[ia] = e.w;
            adj[ib] = a; aw[ib] = e.w;
        }

        order_.reserve(N);
        std::vector<char> seen(N, 0);
        order_.push_back(0);
        seen[0] = 1;
        for (std::size_t qi = 0; qi < order_.size(); ++qi) {
            const std::size_t u = static_cast<std::size_t>(order_[qi]);
            for (int i = head[u]; i < head[u + 1]; ++i) {
                const int w = adj[static_cast<std::size_t>(i)];
                const std::size_t sw = static_cast<std::size_t>(w);
                if (!seen[sw]) {
                    seen[sw] = 1;
                    parent_[sw] = static_cast<int>(u);
                    wpar_[sw] = aw[static_cast<std::size_t>(i)];
                    order_.push_back(w);
                }
            }
        }
    }

    void compute() {
        const std::size_t N = static_cast<std::size_t>(n_);
        f0_.assign(N, {});
        f1_.assign(N, {});
        for (std::size_t v = 0; v < N; ++v) {
            f0_[v] = {0};    // 0 edges, unmatched -- weight 0.
            f1_[v] = {NEG};  // matched needs >= 1 edge -- column 0 unreachable.
        }

        // fold every non-root vertex into its parent, deepest first.
        for (std::size_t idx = order_.size(); idx-- > 1;) {
            const std::size_t u = static_cast<std::size_t>(order_[idx]);
            const std::size_t p = static_cast<std::size_t>(parent_[u]);
            fold(p, u, wpar_[u]);
        }
    }

    // fold child c (final) into parent p across edge weight w.
    void fold(std::size_t p, std::size_t c, i64 w) {
        const std::vector<i64>& p0 = f0_[p];
        const std::vector<i64>& p1 = f1_[p];
        const std::vector<i64>& c0 = f0_[c];
        const std::vector<i64>& c1 = f1_[c];
        const std::size_t Lp = p0.size();
        const std::size_t Lc = c0.size();

        std::vector<i64> n0(Lp + Lc, NEG);  // widen by the child's span (+ the new edge).
        std::vector<i64> n1(Lp + Lc, NEG);

        for (std::size_t a = 0; a < Lp; ++a) {
            const i64 pa0 = p0[a];
            const i64 pa1 = p1[a];
            if (pa0 <= NEG && pa1 <= NEG) continue;
            for (std::size_t b = 0; b < Lc; ++b) {
                // c's contribution when edge (p,c) is skipped: c's own best.
                const i64 cb = std::max(c0[b], c1[b]);
                if (cb > NEG) {
                    if (pa0 > NEG) relax(n0[a + b], pa0 + cb);       // p stays free.
                    if (pa1 > NEG) relax(n1[a + b], pa1 + cb);       // p already matched.
                }
                // take edge (p,c): p free, c free, one more edge, +w.
                if (pa0 > NEG && c0[b] > NEG)
                    relax(n1[a + b + 1], pa0 + c0[b] + w);
            }
        }

        f0_[p].swap(n0);
        f1_[p].swap(n1);
    }

    static void relax(i64& slot, i64 cand) { if (cand > slot) slot = cand; }

    int n_;
    int root_ = 0;
    std::vector<int> order_;
    std::vector<int> parent_;
    std::vector<i64> wpar_;
    std::vector<std::vector<i64>> f0_;
    std::vector<std::vector<i64>> f1_;
};

// convenience: exact g(K), or -1.
inline i64 reference(int K, int n, const std::vector<Edge>& edges) {
    return TreeMatchRef(n, edges).at(K);
}

}  // namespace p0033
