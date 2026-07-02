#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0009 {

// segment tree beats, the light version. one array, three ops:
//   range_sum(l, r)   -- sum of a[l..r]
//   range_mod(l, r, x) -- a[i] %= x for every i in [l, r]
//   point_set(k, x)   -- a[k] = x
//
// each node carries two aggregates over its segment: the sum, and the max. sum
// answers the query; max is what makes the mod cheap.
//
// THE PRUNE. a[i] %= x changes a[i] only when a[i] >= x. so if the largest value
// in a whole subtree is already < x, the mod touches nothing there -- we return
// at the node and never walk in. that one comparison collapses most of the work:
// a range-mod only descends into the parts that actually still hold a value big
// enough to reduce.
//
// WHY IT'S FAST -- the potential argument (amortized, not per-op). when a[i] does
// get modded, a[i] >= x, and a[i] % x < x <= a[i]. sharper: a % x < a/2. proof --
// if x <= a/2 then a%x < x <= a/2; if x > a/2 then a >= x means a%x = a - x < a/2.
// so every real mod at least halves the value. a value starting <= V can be
// halved only O(log V) times before it hits 0 and stops changing. point_set is
// the only thing that lifts a value back up, and it lifts one element by at most
// V, funding another O(log V) halvings. so across the whole run the number of
// leaf-level mods that actually fire is O((n + q_set) log V), q_set = # of sets.
// each such firing costs the O(log n) descent that reached it. every range-mod
// also pays an unavoidable O(log n) just to test the max on its canonical nodes.
// net: O((n + m log n) log V) amortized over m ops, V the value ceiling (1e9).
//
// range_sum and point_set are the ordinary segment-tree O(log n) each, worst
// case -- no amortization hiding in those two.
//
// LOWER BOUND, honest. the sum/update half of this is dynamic partial sums:
// point updates plus prefix/range sums. Patrascu-Demaine put that at
// Omega(log n) per op amortized in the cell-probe model -- you cannot beat
// O(log n) for the query/update mix on arbitrary ranges, so Omega(m log n)
// total is forced regardless of the mod trick. the mod op adds only the
// amortized log V factor on top, and it is amortized: a single range-mod can, in
// isolation, descend to Theta(n) leaves right after a set repopulates them. the
// bound is a run-total, charged against the potential, never a per-op guarantee.
class Beats {
public:
    // build the tree over a[0..n-1]. values are 0-based internally; the driver
    // converts the statement's 1-based indices at the boundary.
    void build(const std::vector<std::int64_t>& a) {
        n_ = static_cast<int>(a.size());
        std::size_t cap = static_cast<std::size_t>(4 * (n_ > 0 ? n_ : 1));
        sum_.assign(cap, 0);
        max_.assign(cap, 0);
        if (n_ > 0) build(1, 0, n_ - 1, a);
    }

    // a[k] = x. an ordinary point update -- O(log n).
    void point_set(int k, std::int64_t x) { point_set(1, 0, n_ - 1, k, x); }

    // a[i] %= x for i in [l, r], inclusive. amortized cheap via the max prune.
    void range_mod(int l, int r, std::int64_t x) {
        range_mod(1, 0, n_ - 1, l, r, x);
    }

    // sum of a[l..r], inclusive. O(log n). caller keeps this in int64 -- a full
    // array of 1e9 over 1e5 cells is 1e14, well past 32 bits.
    std::int64_t range_sum(int l, int r) {
        return range_sum(1, 0, n_ - 1, l, r);
    }

private:
    int n_ = 0;
    std::vector<std::int64_t> sum_;  // segment sums, node-indexed.
    std::vector<std::int64_t> max_;  // segment maxima, node-indexed.

    // recombine a node from its two children.
    void pull(int node) {
        int lc = 2 * node, rc = lc + 1;
        sum_[node] = sum_[lc] + sum_[rc];
        max_[node] = std::max(max_[lc], max_[rc]);
    }

    void build(int node, int lo, int hi, const std::vector<std::int64_t>& a) {
        if (lo == hi) {
            sum_[node] = max_[node] = a[static_cast<std::size_t>(lo)];
            return;
        }
        int mid = (lo + hi) / 2;
        build(2 * node, lo, mid, a);
        build(2 * node + 1, mid + 1, hi, a);
        pull(node);
    }

    void point_set(int node, int lo, int hi, int k, std::int64_t x) {
        if (lo == hi) {  // reached the leaf -- overwrite both aggregates.
            sum_[node] = max_[node] = x;
            return;
        }
        int mid = (lo + hi) / 2;
        if (k <= mid)
            point_set(2 * node, lo, mid, k, x);
        else
            point_set(2 * node + 1, mid + 1, hi, k, x);
        pull(node);
    }

    void range_mod(int node, int lo, int hi, int l, int r, std::int64_t x) {
        if (r < lo || hi < l) return;   // segment sits outside [l, r] -- skip.
        if (max_[node] < x) return;     // THE PRUNE: every value here is already
                                        // < x, so %= x is a no-op on the whole
                                        // subtree, in-range part included.
        if (lo == hi) {                 // a single element that is >= x: reduce.
            sum_[node] %= x;
            max_[node] = sum_[node];
            return;
        }
        int mid = (lo + hi) / 2;
        range_mod(2 * node, lo, mid, l, r, x);
        range_mod(2 * node + 1, mid + 1, hi, l, r, x);
        pull(node);
    }

    std::int64_t range_sum(int node, int lo, int hi, int l, int r) {
        if (r < lo || hi < l) return 0;             // disjoint.
        if (l <= lo && hi <= r) return sum_[node];  // fully covered.
        int mid = (lo + hi) / 2;
        return range_sum(2 * node, lo, mid, l, r) +
               range_sum(2 * node + 1, mid + 1, hi, l, r);
    }
};

}  // namespace p0009
