#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0033 {

using i64 = std::int64_t;
using i128 = __int128;  // penalized values need the headroom -- see below.

// one weighted tree edge, 1-indexed endpoints as they arrive on stdin.
struct Edge {
    int u;
    int v;
    i64 w;
};

// exactly-K maximum-weight matching on a tree, by the ALIENS TRICK.
//
// a matching is a set of edges sharing no endpoint. g(t) = max weight of a
// matching using EXACTLY t edges. we want g(K), or -1 when a t=K matching can't
// exist (K past the tree's maximum matching cardinality). K=0 is 0.
//
// the naive exact route is the tree-knapsack dp f[v][t][matched] -- O(n*K)
// states, and at K ~ n/2 that's ~n^2/2 cells, 2e10 for n=2e5. dead. the aliens
// trick charges a fixed toll lambda per matched edge, solves the toll-adjusted
// problem in ONE linear tree pass, and binary-searches lambda to land on exactly
// K edges. O(n log(n*maxW)). that's the whole game.
//
// LOWER BOUND, stated plainly. you must read the n-1 edges -- Omega(n). you
// can't beat that floor. the trick sits at O(n log(n*maxW)), a log over the
// floor, and it buries the O(n*K) knapsack the moment K is large. the log is the
// price of not knowing the right toll up front. the range is n*maxW, not maxW,
// because past the peak of g the marginal edge can COST up to ~(n/2)*maxW, so the
// toll has to reach that far negative -- hence log(n*maxW) = log n + log maxW.
//
// ===================================================================
//  CONCAVITY -- the fact the whole trick rests on
// ===================================================================
//
// claim: g is concave. 2*g(t) >= g(t+1) + g(t-1) for every t. the marginal
// g(t)-g(t-1) is non-increasing -- the t-th edge you're forced to add is worth
// no more than the (t-1)-th. diminishing returns, provable by exchange.
//
// proof. take an optimal (t+1)-matching A and an optimal (t-1)-matching B. i
// build two matchings, each of size t, whose combined weight is exactly
// w(A)+w(B). that forces 2*g(t) >= w(A)+w(B) = g(t+1)+g(t-1).
//
// let C = A ∩ B -- edges common to both. write A = C ∪ A', B = C ∪ B'. the edge
// set A' ∪ B' touches only vertices C doesn't, and every vertex has degree <= 2
// in it (<=1 from A', <=1 from B'), alternating A'/B'. in a tree there are no
// cycles, so A' ∪ B' is a disjoint union of ALTERNATING PATHS. count edges:
//     |A'| = (t+1) - |C|,   |B'| = (t-1) - |C|,   |A'|+|B'| = 2t - 2|C|.
//
// 2-color the edges of those paths properly -- adjacent edges get different
// colors. each color class is a matching (no shared endpoints) and avoids C's
// vertices, so color-class ∪ C is a legal matching. the two classes together use
// every edge of A' ∪ B' exactly once, total 2t-2|C| edges, so they split as
// (t-|C|) and (t-|C|) -- the surplus paths (two more A'-edges overall) let us
// assign the odd path-ends so both classes come out equal. add C back to each:
// two matchings of size t. weight is CONSERVED -- C counts in both (2*w(C)) and
// A'∪B' is used once (w(A')+w(B')), summing to w(A)+w(B). done.
//
// concavity means the slopes s_t = g(t)-g(t-1) are non-increasing: g rises,
// peaks, then falls. UNLIKE "at most K", exactly-K must report g(K) even PAST the
// peak, where slopes go negative -- that's why the toll here also goes negative.
//
// ===================================================================
//  THE LAGRANGIAN -- one linear tree pass per toll
// ===================================================================
//
// charge a toll lambda per matched edge. solve the unconstrained
//     V(lambda) = max over t of ( g(t) - lambda*t )
// with a rooted tree dp, each matched edge worth (w_e - lambda):
//     dp[v][0] = best penalized value for subtree(v), v NOT matched to a child.
//     dp[v][1] = best penalized value for subtree(v), v matched to one child.
// fold children in one at a time (subtree knapsack, but tracking a single toll-
// adjusted scalar instead of a whole array). the root's max(dp[0],dp[1]) = V.
//
// because g is concave, g(t)-lambda*t is concave in t, maximized on the interval
// { t : s_t >= lambda >= s_{t+1} }. raise lambda, that interval slides left --
// the toll tunes the edge count. the conjugate identity (concave g, integer
// breakpoints):
//     g(K) = min over lambda of ( V(lambda) + lambda*K ).
// weak side is free: V(lambda) >= g(K) - lambda*K for every lambda. equality
// needs a lambda whose optimal-t interval CONTAINS K -- a subgradient at K.
//
// ===================================================================
//  PLATEAU / TIE HANDLING -- the bug-prone spot
// ===================================================================
//
// slopes s_t are integers (edge weights are), non-increasing. at integer lambda
// the optimal t is the interval { t : s_{t+1} <= lambda <= s_t }. when lambda
// equals some slope, several t tie at the same penalized value -- a plateau --
// and the edge count JUMPS. a single lambda can be optimal for a whole run of t
// that straddles K. we must still read off exactly K.
//
// fix: the penalized dp breaks value ties by MINIMIZING edge count. that minimum
// is t_lo(lambda) = |{ t : s_t > lambda }|, the low end of the plateau, and it's
// non-increasing in lambda. binary-search the SMALLEST lambda L with
// t_lo(L) <= K. then:
//     t_lo(L)   = |{ s_t > L   }| <= K            (L satisfies the search)
//     t_lo(L-1) = |{ s_t > L-1 }| = |{ s_t >= L }| >= K+1
// so |{s_t > L}| <= K < |{s_t >= L}| -- the slopes equal to L bracket K. that
// gives s_{K+1} <= L <= s_K: L is a valid subgradient at K even when the plateau
// jumped clean over K. equality holds, g(K) = V(L) + L*K. exactly the array-
// aliens machinery, now over a tree.
//
// FEASIBILITY / the -1 boundary. push lambda to its most negative -- every edge
// becomes pure reward, so the dp takes a MAXIMUM-cardinality matching, and its
// count is maxCard, the largest matching the tree admits. if K > maxCard no
// exactly-K matching exists -- return -1. K=0 returns 0 up front. everything in
// between rides the binary search above.
//
// INT64 WATCH. g(K) fits int64 (K*maxW <= 1e14). but a single penalized value is
// (sum of edge weights) - lambda*(count), and past the peak lambda reaches
// -~(n/2)*maxW while count reaches ~n/2 -- their product blows past 2^63. so the
// dp carries values in __int128. only the final g(K), which fits, comes back as
// int64.

struct Penalized {
    i128 value;  // V(lambda) = max_t ( g(t) - lambda*t )
    i64 count;   // t_lo(lambda) = fewest edges reaching that value
};

// builds the tree once (CSR + a parent-before-child order + each vertex's edge
// weight to its parent), then answers penalized(lambda) in O(n) with no realloc.
// solve() drives the toll binary search. same shape as 0019's TreePaths.
class TreeMatch {
public:
    TreeMatch(int n, const std::vector<Edge>& edges) : n_(n < 0 ? 0 : n) {
        build(edges);
    }

    // one toll-adjusted sweep. returns V(lambda) and the MINIMUM edge count that
    // reaches it -- value ties break toward fewer edges, so that count is exactly
    // t_lo(lambda), the plateau floor the binary search steers on. O(n).
    Penalized penalized(i128 lambda) const {
        const std::size_t N = static_cast<std::size_t>(n_);
        // reset every vertex to its base state: unmatched with nothing chosen,
        // matched unreachable until an edge creates it.
        for (std::size_t v = 0; v < N; ++v) {
            val0_[v] = 0;   cnt0_[v] = 0;
            val1_[v] = NEG; cnt1_[v] = 0;
        }

        // deepest first: children come AFTER their parent in the BFS order_, so
        // reverse order visits every vertex only once all its own children are
        // final. we then FOLD the finished vertex into its parent. no recursion --
        // a bamboo of 2e5 vertices would blow the call stack.
        for (std::size_t idx = order_.size(); idx-- > 1;) {
            const std::size_t u = static_cast<std::size_t>(order_[idx]);
            const std::size_t p = static_cast<std::size_t>(parent_[u]);

            // childBest: u's contribution when the edge (p,u) is NOT taken -- the
            // better of u matched / unmatched, ties toward fewer edges. always
            // reachable (val0_ starts at 0 and only grows).
            i128 cbV;
            i64 cbC;
            if (val1_[u] > val0_[u] ||
                (val1_[u] == val0_[u] && cnt1_[u] < cnt0_[u])) {
                cbV = val1_[u]; cbC = cnt1_[u];
            } else {
                cbV = val0_[u]; cbC = cnt0_[u];
            }

            // fold into p using p's state BEFORE this child (subtree knapsack).
            const i128 o0V = val0_[p]; const i64 o0C = cnt0_[p];
            const i128 o1V = val1_[p]; const i64 o1C = cnt1_[p];

            // p stays unmatched: absorb u's best.
            val0_[p] = o0V + cbV; cnt0_[p] = o0C + cbC;

            // p matched: either it was already matched to an earlier child and u
            // just rides along (a), or we match p to u now -- p was free, u free
            // (its dp[0]), pay (w - lambda) (b).
            const i128 aV = o1V + cbV;                          const i64 aC = o1C + cbC;
            const i128 bV = o0V + val0_[u] + (static_cast<i128>(wpar_[u]) - lambda);
            const i64  bC = o0C + cnt0_[u] + 1;
            if (aV > bV || (aV == bV && aC <= bC)) {
                val1_[p] = aV; cnt1_[p] = aC;
            } else {
                val1_[p] = bV; cnt1_[p] = bC;
            }
        }

        const std::size_t rt = static_cast<std::size_t>(order_.empty() ? 0 : order_[0]);
        if (N == 0) return {0, 0};
        if (val1_[rt] > val0_[rt] ||
            (val1_[rt] == val0_[rt] && cnt1_[rt] < cnt0_[rt]))
            return {val1_[rt], cnt1_[rt]};
        return {val0_[rt], cnt0_[rt]};
    }

    // max matching cardinality -- the count at the most negative toll, where
    // every edge is pure reward. this is the -1 boundary: K above it is infeasible.
    i64 max_cardinality() const { return penalized(lambda_lo()).count; }

    // g(K): max weight of a matching using EXACTLY K edges, or -1 if impossible.
    i64 solve(int K) const {
        if (K <= 0) return 0;                     // empty matching, weight 0.
        if (n_ <= 1) return -1;                   // no edges -- only K=0 is legal.

        const i64 maxCard = max_cardinality();
        if (static_cast<i64>(K) > maxCard) return -1;  // K past the biggest matching.

        // smallest toll L with t_lo(L) <= K. the predicate is monotone (t_lo is
        // non-increasing in lambda), lambda_hi() always satisfies it, so this is a
        // clean lower_bound. for K == maxCard it lands on lambda_lo, which is still
        // a valid subgradient there. L is a subgradient at K -- plateau argument.
        i128 lo = lambda_lo(), hi = lambda_hi();
        while (lo < hi) {
            const i128 mid = lo + (hi - lo) / 2;   // hi-lo >= 0, floor -- safe negative lo.
            if (penalized(mid).count <= static_cast<i64>(K)) hi = mid;
            else lo = mid + 1;
        }
        const Penalized best = penalized(lo);
        // g(K) = V(L) + L*K. fits int64 by the answer bound; the int128 was only
        // for the intermediate penalized value.
        return static_cast<i64>(best.value + lo * static_cast<i128>(K));
    }

private:
    // toll range. hi: s_1 <= maxW, so t_lo(maxW) = 0 -- covers every K. lo: below
    // every slope, so t_lo = maxCard -- the most negative a marginal can be is
    // ~(n/2)*maxW (one augmenting swap along a path), so -(maxW*n)-1 clears it.
    i128 lambda_hi() const { return static_cast<i128>(eff_w()); }
    i128 lambda_lo() const {
        return -(static_cast<i128>(eff_w()) * static_cast<i128>(n_)) - 1;
    }
    i64 eff_w() const { return maxW_ > 0 ? maxW_ : 1; }  // guard the no-edge tree.

    void build(const std::vector<Edge>& edges) {
        const std::size_t N = static_cast<std::size_t>(n_);
        parent_.assign(N, -1);
        wpar_.assign(N, 0);
        val0_.assign(N, 0); cnt0_.assign(N, 0);
        val1_.assign(N, 0); cnt1_.assign(N, 0);
        order_.clear();
        maxW_ = 0;
        if (n_ == 0) return;

        // CSR adjacency, undirected, weights carried alongside. one build, one
        // allocation, cache-friendly -- same pattern 0019 uses.
        std::vector<int> deg(N, 0);
        for (const Edge& e : edges) {
            ++deg[static_cast<std::size_t>(e.u - 1)];
            ++deg[static_cast<std::size_t>(e.v - 1)];
            if (e.w > maxW_) maxW_ = e.w;
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

        // BFS from root 0: parent_, wpar_, and a parent-before-child order_.
        // iterative on purpose -- never touch the call stack.
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

    static constexpr i128 NEG = -(static_cast<i128>(1) << 100);  // unreachable state.

    int n_;
    i64 maxW_ = 0;
    std::vector<int> order_;             // parent-before-child BFS order.
    std::vector<int> parent_;            // parent_[v] = parent, -1 at the root.
    std::vector<i64> wpar_;              // weight of edge (v, parent_[v]).
    mutable std::vector<i128> val0_;     // dp[v][0] value, scratch for penalized().
    mutable std::vector<i64> cnt0_;      // dp[v][0] min edge count.
    mutable std::vector<i128> val1_;     // dp[v][1] value.
    mutable std::vector<i64> cnt1_;      // dp[v][1] min edge count.
};

// convenience wrapper: build once, answer one K.
inline i64 solve(int K, int n, const std::vector<Edge>& edges) {
    return TreeMatch(n, edges).solve(K);
}

}  // namespace p0033
