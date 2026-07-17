#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0044 {

// range k-th smallest. given a[1..n], answer "l r k" -- the k-th smallest value
// in a[l..r]. the naive answer sorts a[l..r] per query: O(q * n log n). the trick
// that kills the n is a persistent segment tree over the value domain.
//
// THE IDEA -- one segment tree per prefix, all sharing memory.
//   root[i] holds the value-histogram of a[1..i]: a segment tree over the
//   compressed values where leaf v carries "how many of a[1..i] equal value v".
//   root[i] = root[i-1] with a[i]'s value +1 at its leaf. an insert only rewrites
//   the one root-to-leaf path -- O(log n) new nodes -- and points every off-path
//   child straight back at root[i-1]'s node. so all n+1 versions cost O(n log n)
//   nodes total, not O(n^2).
//
// WHY IT ANSWERS A RANGE. root[r] and root[l-1] differ by exactly the inserts for
// positions l..r. subtract them node-for-node and the difference IS the histogram
// of a[l..r] -- no tree stored per range, just two versions and a subtraction.
// descend both roots in parallel: at a node, the count landing in the left half
// over a[l..r] is cnt[left(root[r])] - cnt[left(root[l-1])]. if k <= that, the
// k-th smallest lives left; else drop those and go right with k -= that. the leaf
// you reach is the k-th smallest value.
//
// COST. build is n inserts, O(n log n). each query is one parallel root-to-leaf
// descent, O(log n). total O((n + q) log n) time, O(n log n) memory.
//
// LOWER BOUND, honest. you must read the array and the queries -- Omega(n + q),
// unavoidable. each query names a value that depends on the whole multiset of
// a[l..r], and locating the k-th order statistic in a size-m histogram is a
// descent of Omega(log m) comparisons in the comparison/decision-tree sense; the
// tree meets it at O(log n) per query. the naive sort-per-query is O(q n log n) --
// persistence buys back the factor of n by never rebuilding a[l..r], only
// differencing two prefixes that already exist.
class PersistentKth {
public:
    // build all n+1 versions over a[0..n-1] (0-based storage; the driver hands us
    // the statement's values in order). node 0 is the shared null: cnt 0, both
    // children point at itself, so descending an absent subtree stays absent and
    // the count-subtraction sees zeros with no special case.
    void build(const std::vector<std::int64_t>& a) {
        n_ = static_cast<int>(a.size());

        // coordinate compression. duplicates collapse to one leaf -- a value that
        // shows up three times is three +1s at the same leaf, cnt 3, not three
        // leaves. that keeps the domain at "distinct values" and the histogram
        // honest.
        sorted_.assign(a.begin(), a.end());
        std::sort(sorted_.begin(), sorted_.end());
        sorted_.erase(std::unique(sorted_.begin(), sorted_.end()), sorted_.end());
        m_ = static_cast<int>(sorted_.size());
        if (m_ == 0) m_ = 1;  // n==0 degenerate -- keep a one-leaf domain.

        // exact node budget: the null node, plus one root-to-leaf path per insert.
        // path length = tree height = floor(log2(m-1)) + 1 leaves-and-internals,
        // bounded by bits+1. reserve it so no push_back ever reallocates mid-build.
        int bits = 0;
        while ((1 << bits) < m_) ++bits;
        const std::size_t per_insert = static_cast<std::size_t>(bits) + 2;
        lc_.reserve(static_cast<std::size_t>(n_) * per_insert + 2);
        rc_.reserve(static_cast<std::size_t>(n_) * per_insert + 2);
        cnt_.reserve(static_cast<std::size_t>(n_) * per_insert + 2);

        // the null node, index 0. self-childed so left(null)=right(null)=null.
        lc_.assign(1, 0);
        rc_.assign(1, 0);
        cnt_.assign(1, 0);

        root_.assign(static_cast<std::size_t>(n_) + 1, 0);
        root_[0] = 0;  // empty prefix a[1..0] -- the null tree.
        for (int i = 0; i < n_; ++i) {
            int pos = compress(a[static_cast<std::size_t>(i)]);
            root_[static_cast<std::size_t>(i) + 1] =
                insert(root_[static_cast<std::size_t>(i)], 0, m_ - 1, pos);
        }
    }

    // k-th smallest value in a[l..r], every index 1-based, 1 <= k <= r-l+1.
    // descend root[r] against root[l-1] -- the difference is a[l..r]'s histogram.
    std::int64_t query(int l, int r, int k) const {
        int vl = root_[static_cast<std::size_t>(l - 1)];
        int vr = root_[static_cast<std::size_t>(r)];
        int lo = 0, hi = m_ - 1;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            // count in the left half over a[l..r] = left(root[r]) - left(root[l-1]).
            int left_cnt = cnt_[static_cast<std::size_t>(lc_[static_cast<std::size_t>(vr)])] -
                           cnt_[static_cast<std::size_t>(lc_[static_cast<std::size_t>(vl)])];
            if (k <= left_cnt) {  // the k-th smallest sits in the left half.
                vl = lc_[static_cast<std::size_t>(vl)];
                vr = lc_[static_cast<std::size_t>(vr)];
                hi = mid;
            } else {  // skip the whole left half, hunt the (k - left_cnt)-th right.
                k -= left_cnt;
                vl = rc_[static_cast<std::size_t>(vl)];
                vr = rc_[static_cast<std::size_t>(vr)];
                lo = mid + 1;
            }
        }
        return sorted_[static_cast<std::size_t>(lo)];  // leaf reached -- decompress.
    }

private:
    int n_ = 0;
    int m_ = 0;                        // distinct-value domain size.
    std::vector<std::int64_t> sorted_; // compressed value -> original value.
    std::vector<int> root_;            // root_[i] = version after a[1..i].
    std::vector<int> lc_, rc_;         // persistent node children, index 0 = null.
    std::vector<int> cnt_;             // subtree count; n <= 1e5 fits an int.

    // value -> compressed index. lower_bound is exact -- every queried value came
    // from the array, so it's always present.
    int compress(std::int64_t v) const {
        return static_cast<int>(
            std::lower_bound(sorted_.begin(), sorted_.end(), v) - sorted_.begin());
    }

    // clone the path from prev, +1 at leaf pos. one new node per level; the child
    // we don't walk into is aliased straight to prev's -- that aliasing is the
    // whole memory win. returns the new node's index.
    int insert(int prev, int lo, int hi, int pos) {
        int cur = alloc();
        lc_[static_cast<std::size_t>(cur)] = lc_[static_cast<std::size_t>(prev)];
        rc_[static_cast<std::size_t>(cur)] = rc_[static_cast<std::size_t>(prev)];
        cnt_[static_cast<std::size_t>(cur)] = cnt_[static_cast<std::size_t>(prev)] + 1;
        if (lo == hi) return cur;  // leaf -- the +1 already landed in cnt.
        int mid = lo + (hi - lo) / 2;
        if (pos <= mid) {
            lc_[static_cast<std::size_t>(cur)] =
                insert(lc_[static_cast<std::size_t>(prev)], lo, mid, pos);
        } else {
            rc_[static_cast<std::size_t>(cur)] =
                insert(rc_[static_cast<std::size_t>(prev)], mid + 1, hi, pos);
        }
        return cur;
    }

    // hand out the next node index. capacity was reserved in build, so this never
    // reallocates -- and even if it did, every reference here is an index, not a
    // pointer, so a move would be harmless.
    int alloc() {
        int id = static_cast<int>(lc_.size());
        lc_.push_back(0);
        rc_.push_back(0);
        cnt_.push_back(0);
        return id;
    }
};

}  // namespace p0044
