#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace p0019 {

// codeforces 1039D -- you are given a tree. n vertices, unrooted tree.
// a "path of length k" is exactly k vertices forming a simple path. for every
// k from 1 to n, report f(k) = the maximum number of vertex-disjoint length-k
// paths you can lay on the tree at once.
//
// honest lower bound: you must read the n-1 edges (Omega(n)) and emit n answers
// (Omega(n)). so Omega(n) is the floor. the naive route -- run the O(n) greedy
// once per k -- is O(n^2). the sqrt-on-the-answer trick below cuts the number of
// greedy runs from n to O(sqrt(n) * log n), landing at O(n * sqrt(n) * log n).
// the log is the binary search inside each constant run; drop it in the headline
// and people call this "O(n sqrt n)", but the run time carries the log, so that
// is what we state.
//
// ===================================================================
//  part 1 -- the O(n) greedy for a FIXED k, and why it is optimal
// ===================================================================
//
// root anywhere (vertex 0). process vertices children-before-parent. each vertex
// hands its parent a "tail": the length, in vertices, of the longest chain of
// still-unused vertices that runs straight up to it and could be finished into a
// path higher up. a used vertex hands up 0.
//
// at vertex v let d1 >= d2 be the two largest tails among v's children (0 when v
// has fewer than two children). a path that peaks at v goes down one child, up
// through v, down another child, so its vertex count is (x) + 1 + (y) with
// 0 <= x <= d1 and 0 <= y <= d2. it can hit exactly k iff k-1 <= d1 + d2, i.e.
//        1 + d1 + d2 >= k.
// (a single-arm path uses d2 = 0, already covered.) so:
//   * if 1 + d1 + d2 >= k  -> close a path here (answer += 1), hand up tail 0.
//   * otherwise            -> hand up tail 1 + d1 (v extends the longest child).
//
// once v is spent inside a closed path, no chain can pass THROUGH v upward -- the
// only route up is v itself -- so every other child tail below v is wasted. that
// is why "close" always hands up 0, and why the only real decision is close vs.
// don't-close. using the top two tails for the feasibility test is right because
// if the two longest can't reach k, no pair can.
//
// ---- optimality (exchange / induction), the subtle part ----
//
// for a vertex v define greedy's pair (cnt_g(v), tail_g(v)) = (paths closed
// inside v's subtree, tail handed up). claim, by induction over the subtree:
//
//   for EVERY legal packing S of length-k paths lying inside subtree(v), letting
//   cnt(S) be its path count and tail(S) the longest unused straight-down chain
//   it leaves at v (tail 0 if v is used),
//        cnt(S) <= cnt_g(v),  and  cnt(S) == cnt_g(v)  =>  tail(S) <= tail_g(v).
//
// so greedy maximizes closed paths and, among all maximizers, the tail it exposes
// -- the strongest possible offer to the vertex above. at the root the tail is
// irrelevant and cnt_g(root) >= cnt(S) for all S: greedy is optimal.
//
// base (leaf): only packings are {v unused: cnt 0, tail 1} and, when k == 1,
// {v is its own path: cnt 1, tail 0}. greedy closes iff 1 >= k iff k == 1 -- both
// arms match.
//
// step (v with children): S restricts to each child subtree as (cnt_i, tail_i),
// and by IH cnt_i <= cnt_g(c_i) with tail_i <= tail_g(c_i) whenever cnt_i is
// maximal. write C* = sum_i cnt_g(c_i); then sum_i cnt_i <= C*. two cases for v:
//
//   (A) S closes a path peaking at v. it eats two child tails tail_a, tail_b with
//       1 + tail_a + tail_b >= k, and cnt(S) = sum_i cnt_i + 1, tail(S) = 0.
//       - if greedy also closes: cnt_g(v) = C* + 1 >= sum_i cnt_i + 1 = cnt(S).
//       - if greedy does NOT close (1 + t1 + t2 < k, t1,t2 the top greedy tails):
//         should S reach cnt(S) = C* + 1 it would need sum_i cnt_i = C*, i.e.
//         every child maximal, forcing tail_a,tail_b <= their tail_g <= t1,t2 by
//         choice of top two -- then 1 + tail_a + tail_b <= 1 + t1 + t2 < k, so S
//         cannot close at v. contradiction. hence sum_i cnt_i <= C* - 1 and
//         cnt(S) <= C* = cnt_g(v). either way cnt(S) <= cnt_g(v); tail(S)=0 ok.
//
//   (B) v is unused. cnt(S) = sum_i cnt_i <= C* <= cnt_g(v). equality forces
//       greedy to not close (else cnt_g(v) = C*+1 > cnt(S)) and every child
//       maximal, so tail(S) = 1 + max_i tail_i <= 1 + t1 = tail_g(v).
//
// the induction closes. greedy is exactly optimal. the differential test against
// the independent brute oracle hammers this on thousands of trees across all k.
//
// ===================================================================
//  part 2 -- f takes O(sqrt n) distinct values, so O(sqrt n) runs
// ===================================================================
//
// f is NON-INCREASING: from m disjoint (k+1)-vertex paths, delete one endpoint of
// each -- still a simple path, now k vertices, still disjoint. so f(k) >= f(k+1).
//
// f(k) <= n / k: the paths are disjoint and each uses k vertices, so k*f(k) <= n.
//
// count distinct values. for k > sqrt(n): f(k) <= n/k < sqrt(n), so f lands in
// {0,...,floor(sqrt n)} -- at most sqrt(n)+1 distinct values there. for
// k <= sqrt(n): at most sqrt(n) values of k, hence <= sqrt(n) values of f. total
// distinct values of f <= 2*sqrt(n) + 1 = O(sqrt n). since f is non-increasing,
// [1, n] splits into O(sqrt n) maximal runs on which f is constant.
//
// so: walk k left to right. at the start of a run compute val = f(k) with one
// greedy pass, then binary-search the largest r with f(r) == val (== f(r) >= val,
// as f can't exceed val past k), fill [k, r] with val, jump to r+1. O(sqrt n)
// runs, each one binary search of O(log n) greedy passes, each pass O(n).

// holds the tree once (CSR + a parent-before-child order) and answers greedy(k)
// for any k without rebuilding. greedy(k) is O(n); solve() drives the sqrt loop.
class TreePaths {
public:
    TreePaths(int n, const std::vector<std::pair<int, int>>& edges)
        : n_(n < 0 ? 0 : n) {
        build(edges);
    }

    // maximum number of vertex-disjoint length-k paths. O(n), no allocation.
    // best1_/best2_ hold each vertex's top two child tails; we clear a vertex's
    // slots right after reading them, so they are back to 0 for the next call --
    // no per-call fill, the loop pays for itself.
    int greedy(int k) const {
        int count = 0;
        // reverse of a parent-before-child order = children before parents.
        for (std::size_t idx = order_.size(); idx-- > 0;) {
            const int u = order_[idx];
            const std::size_t su = static_cast<std::size_t>(u);
            const int d1 = best1_[su];
            const int d2 = best2_[su];

            int tail;
            if (1 + d1 + d2 >= k) {
                ++count;      // a length-k path peaks here.
                tail = 0;     // v is spent -- nothing passes up through it.
            } else {
                tail = 1 + d1;  // v extends the longest child chain upward.
            }

            best1_[su] = 0;  // reset for the next greedy() call.
            best2_[su] = 0;

            const int p = parent_[su];
            if (p >= 0) {
                const std::size_t sp = static_cast<std::size_t>(p);
                if (tail > best1_[sp]) {
                    best2_[sp] = best1_[sp];
                    best1_[sp] = tail;
                } else if (tail > best2_[sp]) {
                    best2_[sp] = tail;
                }
            }
        }
        return count;
    }

    // f(k) for every k in [1, n], returned as ans[k-1]. drives the sqrt-on-answer
    // loop: one greedy per run start, one binary search per run.
    std::vector<int> solve() const {
        std::vector<int> ans(static_cast<std::size_t>(n_), 0);
        if (n_ == 0) return ans;

        int k = 1;
        while (k <= n_) {
            const int val = greedy(k);
            if (val == 0) {
                // f only drops from here -- every larger k is 0 already. done.
                break;
            }
            // largest r with f(r) == val. f is non-increasing and f(k) == val, so
            // for r >= k the tests f(r) == val and f(r) >= val coincide. and since
            // f(r) == val forces r*val <= n (each of val disjoint paths eats r
            // vertices), the run can't reach past n/val -- clamp hi there. this is
            // what actually earns the sqrt bound: high-value runs sit at small k
            // with a tiny r-window, so their binary search is a call or two, not
            // log n wasted passes over the whole tail.
            int lo = k, hi = n_ / val, r = k;
            while (lo <= hi) {
                const int mid = lo + (hi - lo) / 2;
                if (greedy(mid) >= val) {
                    r = mid;
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
            for (int j = k; j <= r; ++j)
                ans[static_cast<std::size_t>(j - 1)] = val;
            k = r + 1;
        }
        return ans;
    }

private:
    void build(const std::vector<std::pair<int, int>>& edges) {
        best1_.assign(static_cast<std::size_t>(n_), 0);
        best2_.assign(static_cast<std::size_t>(n_), 0);
        parent_.assign(static_cast<std::size_t>(n_), -1);
        if (n_ == 0) return;

        const std::size_t N = static_cast<std::size_t>(n_);

        // ---- CSR adjacency, undirected. one contiguous neighbor array, one
        // allocation, cache-friendly -- the same build 0015 uses. O(n). ----
        std::vector<int> deg(N, 0);
        for (const auto& e : edges) {
            ++deg[static_cast<std::size_t>(e.first - 1)];
            ++deg[static_cast<std::size_t>(e.second - 1)];
        }
        std::vector<int> head(N + 1, 0);
        for (std::size_t v = 0; v < N; ++v) head[v + 1] = head[v] + deg[v];
        std::vector<int> adj(static_cast<std::size_t>(head[N]));
        std::vector<int> cur = head;
        for (const auto& e : edges) {
            const int a = e.first - 1, b = e.second - 1;
            adj[static_cast<std::size_t>(cur[static_cast<std::size_t>(a)]++)] = b;
            adj[static_cast<std::size_t>(cur[static_cast<std::size_t>(b)]++)] = a;
        }

        // ---- BFS from root 0: parent[] and an order with parent before child.
        // iterative on purpose -- a bamboo of 1e5 vertices would blow a recursive
        // DFS stack, so we never touch the call stack. ----
        order_.clear();
        order_.reserve(N);
        std::vector<char> seen(N, 0);
        order_.push_back(0);
        seen[0] = 1;
        for (std::size_t qi = 0; qi < order_.size(); ++qi) {
            const int u = order_[qi];
            for (int i = head[static_cast<std::size_t>(u)];
                 i < head[static_cast<std::size_t>(u) + 1]; ++i) {
                const int w = adj[static_cast<std::size_t>(i)];
                if (!seen[static_cast<std::size_t>(w)]) {
                    seen[static_cast<std::size_t>(w)] = 1;
                    parent_[static_cast<std::size_t>(w)] = u;
                    order_.push_back(w);
                }
            }
        }
    }

    int n_;
    std::vector<int> order_;           // parent-before-child BFS order.
    std::vector<int> parent_;          // parent_[v] = parent, -1 at the root.
    mutable std::vector<int> best1_;   // top child tail, scratch for greedy().
    mutable std::vector<int> best2_;   // second child tail.
};

// convenience wrapper: build once, solve. ans[k-1] = f(k) for k in [1, n].
inline std::vector<int> you_are_given_a_tree(
    int n, const std::vector<std::pair<int, int>>& edges) {
    return TreePaths(n, edges).solve();
}

}  // namespace p0019
