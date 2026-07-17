#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace p0045 {

// segment tree beats -- the Ji Driver tree. four ops on a[0..n-1]:
//   chmin(l, r, x) -- a[i] = min(a[i], x) for i in [l, r]
//   add(l, r, x)   -- a[i] += x           for i in [l, r]   (x may be < 0)
//   query_max(l, r) -- max of a[l..r]
//   query_sum(l, r) -- sum of a[l..r]
//
// the naive chmin walks every cell in range -- O(q n). beats does it in
// O(n log n + q log^2 n) amortized, and the whole trick is one comparison.
//
// EACH NODE CARRIES five numbers over its segment:
//   sum -- the sum. answers query_sum, kept exact at all times.
//   mx  -- the maximum value.
//   se  -- the STRICT second max: the largest value strictly below mx, or
//          -inf if every element equals mx. this is the break condition's fuel.
//   cmx -- how many elements equal mx.
//   add -- a pending additive lazy tag for descendants.
// the chmin has no separate tag. reducing mx IS the tag: a child whose mx tops
// its parent's mx owes a clamp, discovered lazily on push_down.
//
// THE BREAK CONDITION -- chmin(x) at a node fully inside [l, r]:
//   x >= mx      -- nothing here exceeds x. no-op, return.
//   se < x < mx  -- only the max elements exceed x. they all drop to x and
//                   nothing else moves: sum -= (mx - x) * cmx; mx = x. return
//                   WITHOUT recursing. this is the "beat" -- the node collapses
//                   a whole band of equal maxima in O(1).
//   x <= se      -- x cuts below the second max, so >= two distinct values must
//                   change. only now do we push down and RECURSE.
// recurse-only-when-x<=se is not an optimization bolted onto a correct tree. it
// IS the tree. drop it and chmin becomes O(range) and the amortized bound dies.
//
// WHY THE BREAK IS AMORTIZED-CHEAP -- the potential argument. define
//   Phi = sum over all nodes of (number of DISTINCT values in that node's band).
// Phi starts at O(n log n): a leaf holds one distinct value, and summed up the
// tree that is O(n) per level over log n levels. now watch each op move it:
//   - a beat (se < x < mx) merges the top band into the second band at that
//     node: one fewer distinct value there. every node that beats pays for
//     itself by dropping Phi by >= 1. the descent that reached it is O(log n).
//   - when chmin recurses (x <= se), the recursion bottoms out either at a beat
//     or at a node the range only partially covers. the partial-cover nodes are
//     the O(log n) canonical boundary nodes -- charged to the op, not to Phi.
//     so a single chmin does O(log n) unavoidable boundary work plus O(log n)
//     per unit of Phi it burns.
//   - add raises Phi by at most O(log n): it can pull apart values that a chmin
//     had merged, but only along the O(log n) canonical nodes of its range, each
//     re-splitting by at most one. so q adds inject O(q log n) potential total.
// total beats over the run <= starting Phi + injected Phi = O(n log n + q log n),
// each beat/recursion step costing O(log n) to reach: O(n log n + q log^2 n).
// the log^2 is the descent (log n) times the number of levels a chmin can spend
// potential across (log n). sum and max queries are plain O(log n), no amort.
//
// LOWER BOUND, honest. reading the input is Omega(n + q) -- unavoidable. the
// sum half alone (point-reachable updates + arbitrary range sums) is dynamic
// partial sums, which Patrascu-Demaine pin at Omega(log n) per op amortized in
// the cell-probe model. so Omega(q log n) is forced no matter how clever the
// chmin is. beats sits one log factor above that floor on the chmin op, and the
// gap is amortized, not per-op: a single chmin can, right after an add re-splits
// a settled array, descend to Theta(n) leaves. the bound is a run-total charged
// to Phi, never a promise about any one call.
//
// OVERFLOW. values live in [-1e9 - q*1e9, 1e9 + q*1e9] ~ [-1e14, 1e14] at
// q = 1e5. a single node's sum over n = 1e5 such cells reaches ~1e19, past the
// int64 ceiling of ~9.2e18. the statement therefore CAPS the reachable
// magnitude: any element stays within |a_i| <= 9e13 across the run (equivalently
// the adds don't all pile onto one cell for the full q), which keeps every
// partial sum under n * 9e13 = 9e18 < 2^63. inside that cap int64 is exact; the
// differential tests stay in a tiny value band so they never approach it.
class Beats {
public:
    // build over a[0..n-1]. 0-based internally; the driver shifts the
    // statement's 1-based indices at the boundary.
    void build(const std::vector<std::int64_t>& a) {
        n_ = static_cast<int>(a.size());
        std::size_t cap = static_cast<std::size_t>(4 * (n_ > 0 ? n_ : 1));
        sum_.assign(cap, 0);
        mx_.assign(cap, 0);
        se_.assign(cap, kNegInf);
        cmx_.assign(cap, 0);
        add_.assign(cap, 0);
        len_.assign(cap, 0);
        if (n_ > 0) build(1, 0, n_ - 1, a);
    }

    // a[i] = min(a[i], x) on [l, r]. amortized cheap via the break condition.
    void chmin(int l, int r, std::int64_t x) { chmin(1, 0, n_ - 1, l, r, x); }

    // a[i] += x on [l, r]. x may be negative. ordinary lazy add.
    void add(int l, int r, std::int64_t x) { add(1, 0, n_ - 1, l, r, x); }

    // max of a[l..r]. O(log n).
    std::int64_t query_max(int l, int r) {
        return query_max(1, 0, n_ - 1, l, r);
    }

    // sum of a[l..r]. O(log n). caller keeps this in int64.
    std::int64_t query_sum(int l, int r) {
        return query_sum(1, 0, n_ - 1, l, r);
    }

private:
    // -inf sentinel for se. quartered so add can shift a REAL se by up to 1e14
    // without ever colliding with it, and the guard below leaves the sentinel
    // itself untouched -- se == kNegInf means "no second value", not a number.
    static constexpr std::int64_t kNegInf =
        std::numeric_limits<std::int64_t>::min() / 4;

    int n_ = 0;
    std::vector<std::int64_t> sum_;   // segment sum, exact.
    std::vector<std::int64_t> mx_;    // segment max.
    std::vector<std::int64_t> se_;    // strict second max, or kNegInf.
    std::vector<std::int64_t> cmx_;   // count of elements equal to mx.
    std::vector<std::int64_t> add_;   // pending additive tag for descendants.
    std::vector<std::int64_t> len_;   // element count in the segment.

    // shift a whole node by v. mx and se ride the shift together so their gap --
    // the thing the break condition tests -- is preserved. se stays the sentinel
    // if it was one: adding to -inf would corrupt "no second value".
    void apply_add(int node, std::int64_t v) {
        sum_[node] += v * len_[node];
        mx_[node] += v;
        if (se_[node] != kNegInf) se_[node] += v;
        add_[node] += v;
    }

    // clamp the max band down to v. PRECONDITION se < v < mx, so only the cmx
    // top elements move and se/cmx are untouched -- that is what makes it O(1).
    void apply_cap(int node, std::int64_t v) {
        sum_[node] -= (mx_[node] - v) * cmx_[node];
        mx_[node] = v;
    }

    // recombine a node from its children. mx is the larger child max; se is the
    // largest value strictly below it, which is the max of both children's se
    // plus whichever child max is itself below the winning mx.
    void pull(int node) {
        int L = 2 * node, R = L + 1;
        sum_[node] = sum_[L] + sum_[R];
        if (mx_[L] == mx_[R]) {
            mx_[node] = mx_[L];
            cmx_[node] = cmx_[L] + cmx_[R];
            se_[node] = std::max(se_[L], se_[R]);
        } else if (mx_[L] > mx_[R]) {
            mx_[node] = mx_[L];
            cmx_[node] = cmx_[L];
            se_[node] = std::max(se_[L], mx_[R]);
        } else {
            mx_[node] = mx_[R];
            cmx_[node] = cmx_[R];
            se_[node] = std::max(se_[R], mx_[L]);
        }
    }

    // hand this node's two pending tags to its children. ADD FIRST, then the
    // cap -- the add and cap were recorded in that time order (a cap only ever
    // reduces the already-added mx), so the children must see them in it. after
    // the add both children sit in the parent's coordinate frame, so a child mx
    // above the parent mx is exactly the clamp the parent's chmin still owes.
    // se_[child] < mx_[parent] is guaranteed: the parent's cap dropped mx below
    // the OLD parent se, and old parent se >= every child se (pull maintains it),
    // and the add shifts both sides equally -- so apply_cap's precondition holds.
    void push_down(int node) {
        int L = 2 * node, R = L + 1;
        if (add_[node] != 0) {
            apply_add(L, add_[node]);
            apply_add(R, add_[node]);
            add_[node] = 0;
        }
        if (mx_[L] > mx_[node]) apply_cap(L, mx_[node]);
        if (mx_[R] > mx_[node]) apply_cap(R, mx_[node]);
    }

    void build(int node, int lo, int hi, const std::vector<std::int64_t>& a) {
        len_[node] = hi - lo + 1;
        if (lo == hi) {
            sum_[node] = mx_[node] = a[static_cast<std::size_t>(lo)];
            se_[node] = kNegInf;  // a single value has no second value.
            cmx_[node] = 1;
            return;
        }
        int mid = (lo + hi) / 2;
        build(2 * node, lo, mid, a);
        build(2 * node + 1, mid + 1, hi, a);
        pull(node);
    }

    void chmin(int node, int lo, int hi, int l, int r, std::int64_t x) {
        // out of range, or nothing here tops x -- either way a no-op. the
        // x >= mx test is safe even on a partially covered node: every element
        // in the whole subtree is <= mx <= x, so min(a[i], x) = a[i] for all.
        if (r < lo || hi < l || x >= mx_[node]) return;
        if (l <= lo && hi <= r && se_[node] < x) {
            apply_cap(node, x);  // THE BEAT: se < x < mx, collapse the top band.
            return;
        }
        // x <= se (or partial cover) -- two distinct values must move, recurse.
        push_down(node);
        int mid = (lo + hi) / 2;
        chmin(2 * node, lo, mid, l, r, x);
        chmin(2 * node + 1, mid + 1, hi, l, r, x);
        pull(node);
    }

    void add(int node, int lo, int hi, int l, int r, std::int64_t x) {
        if (r < lo || hi < l) return;
        if (l <= lo && hi <= r) {
            apply_add(node, x);
            return;
        }
        push_down(node);
        int mid = (lo + hi) / 2;
        add(2 * node, lo, mid, l, r, x);
        add(2 * node + 1, mid + 1, hi, l, r, x);
        pull(node);
    }

    std::int64_t query_max(int node, int lo, int hi, int l, int r) {
        if (r < lo || hi < l) return kNegInf;      // disjoint.
        if (l <= lo && hi <= r) return mx_[node];   // fully covered.
        push_down(node);  // descendants may carry stale tags -- settle them.
        int mid = (lo + hi) / 2;
        return std::max(query_max(2 * node, lo, mid, l, r),
                        query_max(2 * node + 1, mid + 1, hi, l, r));
    }

    std::int64_t query_sum(int node, int lo, int hi, int l, int r) {
        if (r < lo || hi < l) return 0;                 // disjoint.
        if (l <= lo && hi <= r) return sum_[node];      // fully covered.
        push_down(node);  // sum is exact per node, but grandchildren under an
                          // unpushed tag are stale -- push before we read them.
        int mid = (lo + hi) / 2;
        return query_sum(2 * node, lo, mid, l, r) +
               query_sum(2 * node + 1, mid + 1, hi, l, r);
    }
};

}  // namespace p0045
