#pragma once

#include <algorithm>
#include <vector>

namespace p0008 {

// the dumb, obviously-correct oracle -- O(n^2 k). same recurrence as the fast
// one, dp[j][i] = max over the last cut l of dp[j-1][l] + distinct(l+1, i), but
// with zero cleverness in the transition to be wrong in: for each (j, i) it
// walks l from i-1 down to 0 and recomputes distinct(l+1, i) by hand, growing a
// fresh boolean "seen" set one position at a time. no windows, no monotonicity,
// no divide and conquer. only sane for small n -- which is exactly where the
// differential test lives.
//
// types is 0-indexed of length n; the 1-indexed piece a[l+1 .. i] is the
// 0-indexed slice types[l .. i-1], so pushing position l+1 into the window means
// reading types[l].
inline int reference(const std::vector<int>& types, int k) {
    const int n = static_cast<int>(types.size());
    const int kNeg = -1'000'000'000;

    std::vector<int> prev(static_cast<std::size_t>(n) + 1, kNeg);
    std::vector<int> cur(static_cast<std::size_t>(n) + 1, kNeg);
    prev[0] = 0;  // zero pieces, zero cakes, zero value.

    for (int j = 1; j <= k; ++j) {
        std::fill(cur.begin(), cur.end(), kNeg);
        for (int i = 1; i <= n; ++i) {
            // distinct(l+1, i) as l shrinks: start empty at l = i-1 (which adds
            // only position i) and add one cake per step downward. a plain
            // per-type boolean, rebuilt from scratch every i -- can't be subtle.
            std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
            int distinct = 0;
            int best = kNeg;
            for (int l = i - 1; l >= 0; --l) {
                int t = types[static_cast<std::size_t>(l)];  // cake at 1-indexed l+1.
                if (!seen[static_cast<std::size_t>(t)]) {
                    seen[static_cast<std::size_t>(t)] = 1;
                    ++distinct;
                }
                if (prev[static_cast<std::size_t>(l)] != kNeg) {
                    best = std::max(best, prev[static_cast<std::size_t>(l)] + distinct);
                }
            }
            cur[static_cast<std::size_t>(i)] = best;
        }
        prev.swap(cur);
    }
    return prev[static_cast<std::size_t>(n)];
}

}  // namespace p0008
