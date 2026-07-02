#pragma once

#include <algorithm>
#include <vector>

namespace p0008 {

// n cakes in a row, cake i has type a[i]. cut the row into exactly k contiguous
// non-empty pieces. a piece is worth the number of distinct types inside it.
// maximize the sum of the pieces. (codeforces 834D -- n <= 35000, k <= 50.)
//
// the state is honest and unavoidable: dp[j][i] = best total using j pieces to
// cover the first i cakes = max over the last cut l of dp[j-1][l] + d(l+1, i),
// where d(l+1, i) is distinct types in a[l+1 .. i]. that's n*k meaningful cells
// -- every prefix length against every piece count -- so filling the table is
// Omega(nk) no matter how clever the transition is. the naive transition scans
// every l, giving O(n^2 k). the whole game is shaving that inner scan.
//
// two things make it O(n log n) per layer -- O(nk log n) total:
//
//   1. divide and conquer on the argmax. for a fixed j, let opt(i) be the last
//      cut that maximizes dp[j][i]. opt is nondecreasing in i. so compute the
//      middle i first, find its opt, and the left half's answers all sit in
//      [.., opt] while the right half's all sit in [opt, ..]. log n levels, and
//      each level touches O(n) candidate cuts total. that's the log.
//
//   2. distinct(l+1, i) with a movable window. keep counts per type over the
//      current interval [cur_l, cur_r] and a running distinct total. slide an
//      endpoint one step -- add or remove one position -- and the total updates
//      in O(1) (Mo-style). the window persists across the entire layer's
//      recursion, so the endpoints travel O(n log n) steps amortized, not per
//      query. this is the bug-prone part: every add must be paired with the
//      exact remove, and the four while-loops below expand before they shrink so
//      a count never dips below zero.
//
// WHY opt is monotone (the claim the whole method rests on). the cut cost
// w(l, i) = d(l+1, i) obeys the quadrangle inequality: for l1 <= l2 and
// i1 <= i2,
//        w(l1, i1) + w(l2, i2)  >=  w(l1, i2) + w(l2, i1).
// read the gap w(l1, i) - w(l2, i): it counts types living in the extra left
// stretch [l1+1, l2] that do NOT already appear in [l2+1, i] -- the types the
// wider-left window uniquely contributes. grow i and [l2+1, i] swallows more of
// them, so that gap only shrinks. a shrinking gap in i is exactly the
// inequality, and the QI on an additive dp forces the row maxima -- the opt(i)
// -- to march right as i grows. that monotonicity is what lets the divide and
// conquer split the candidate range in half and never miss the true best cut.
//
// values fit int: sum of distinct counts <= sum of piece lengths = n <= 35000.
class Bakery {
public:
    Bakery(const std::vector<int>& types, int k)
        : n_(static_cast<int>(types.size())),
          k_(k),
          a_(static_cast<std::size_t>(n_) + 1),
          cnt_(static_cast<std::size_t>(n_) + 1, 0),
          dp_prev_(static_cast<std::size_t>(n_) + 1, kNeg),
          dp_cur_(static_cast<std::size_t>(n_) + 1, kNeg) {
        // shift to 1-indexed so a_[i] is "the i-th cake" with no off-by-one.
        for (int i = 1; i <= n_; ++i) a_[static_cast<std::size_t>(i)] = types[static_cast<std::size_t>(i - 1)];
    }

    int run() {
        // layer 0: zero pieces cover exactly zero cakes, worth zero. every other
        // prefix length is unreachable with no pieces -- kNeg poisons it.
        std::fill(dp_prev_.begin(), dp_prev_.end(), kNeg);
        dp_prev_[0] = 0;

        for (int j = 1; j <= k_; ++j) {
            std::fill(dp_cur_.begin(), dp_cur_.end(), kNeg);
            // piece count j needs at least j cakes, so i ranges [j, n]; the last
            // cut l ranges [j-1, i-1], hence the seed opt-range [j-1, n-1]. the
            // window carries over from the previous layer -- resetting it would
            // throw away the amortization, and get() moves from wherever it is.
            solve(j, n_, j - 1, n_ - 1);
            dp_prev_.swap(dp_cur_);
        }
        return dp_prev_[static_cast<std::size_t>(n_)];
    }

private:
    // a prefix length can be unreachable (too few cakes for the piece count). a
    // sentinel far below any real value, with headroom so sentinel + n never
    // overflows and never accidentally beats a real candidate.
    static constexpr int kNeg = -1'000'000'000;

    int n_;
    int k_;
    std::vector<int> a_;    // 1-indexed cake types.
    std::vector<int> cnt_;  // cnt_[t] = occurrences of type t in the window.
    std::vector<int> dp_prev_;
    std::vector<int> dp_cur_;

    int distinct_ = 0;  // number of types with cnt_ > 0 -- the window's value.
    int cur_l_ = 1;     // window is a_[cur_l_ .. cur_r_]; starts empty (l > r).
    int cur_r_ = 0;

    // add/remove exactly one boundary position and keep distinct_ exact. a type
    // enters the distinct set the moment its count leaves zero, and leaves the
    // moment its count returns to zero -- nowhere else.
    void add(int pos) {
        if (cnt_[static_cast<std::size_t>(a_[static_cast<std::size_t>(pos)])]++ == 0) ++distinct_;
    }
    void remove(int pos) {
        if (--cnt_[static_cast<std::size_t>(a_[static_cast<std::size_t>(pos)])] == 0) --distinct_;
    }

    // distinct types in a_[l .. r]. expand both ends first, then shrink both --
    // that ordering guarantees we never remove a position the window doesn't
    // currently hold, so no count ever goes negative.
    int window(int l, int r) {
        while (cur_r_ < r) add(++cur_r_);
        while (cur_l_ > l) add(--cur_l_);
        while (cur_r_ > r) remove(cur_r_--);
        while (cur_l_ < l) remove(cur_l_++);
        return distinct_;
    }

    // fill dp_cur_ for prefix lengths [lo, hi], knowing the best last-cut lands
    // in [optl, optr]. solve the midpoint by brute over its slice of the cut
    // range, then recurse -- the found opt splits the cut range for the halves.
    void solve(int lo, int hi, int optl, int optr) {
        if (lo > hi) return;
        int mid = (lo + hi) / 2;

        int best = kNeg;
        int opt = optl;
        // the last cut l sits at most at mid-1 (a piece is non-empty), and at
        // most at optr (monotonicity). dp_prev_[l] finite for all l >= j-1 <=
        // optl, so no sentinel ever wins here.
        int cut_hi = std::min(mid - 1, optr);
        for (int l = optl; l <= cut_hi; ++l) {
            if (dp_prev_[static_cast<std::size_t>(l)] == kNeg) continue;
            int val = dp_prev_[static_cast<std::size_t>(l)] + window(l + 1, mid);
            if (val > best) {
                best = val;
                opt = l;
            }
        }
        dp_cur_[static_cast<std::size_t>(mid)] = best;

        solve(lo, mid - 1, optl, opt);
        solve(mid + 1, hi, opt, optr);
    }
};

// max total value splitting `types` into exactly k contiguous non-empty pieces.
// inputs: types[i] in [1, n]; 1 <= k <= n. output: the maximum, an int in [0, n].
// O(n k log n) time, O(n) space per layer via the two rolling dp rows.
inline int solve(const std::vector<int>& types, int k) {
    if (k <= 0) return 0;  // no pieces, nothing to score -- defensive only.
    Bakery bakery(types, k);
    return bakery.run();
}

}  // namespace p0008
