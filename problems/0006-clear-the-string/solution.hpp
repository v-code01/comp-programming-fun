#pragma once

#include <string>
#include <vector>

namespace p0006 {

// min operations to delete the whole string, where one operation erases a
// contiguous run of one repeated letter and stitches the two sides together.
//
// interval dp. dp[l][r] = cost to clear s[l..r]. look at the leftmost letter
// s[l] and ask how it leaves. two ways, take the cheaper:
//
//   it leaves alone -- one operation just for s[l], then clear the rest:
//     dp[l][r] = dp[l+1][r] + 1
//
//   it leaves glued to a later equal letter s[k] (s[k] == s[l]). before those
//   two can touch, the gap s[l+1..k-1] has to be gone -- that's dp[l+1][k-1].
//   once it's gone s[l] sits flush against s[k], same letter, so it rides out
//   on whatever operation clears s[k..r] for free:
//     dp[l][r] = dp[l+1][k-1] + dp[k][r]
//
// the second form is the one a greedy misses and the one an O(n^2) recurrence
// can't fake -- s[l] may need to reach a letter in the *middle*, not the end.
// "abaa" is the proof: the first 'a' must merge with the trailing "aa" across
// the 'b', costing 2, and any recurrence that only pairs the two endpoints
// reports 3. so we scan every equal k. that scan is the n^3.
//
// base: a single letter is one operation. an empty interval is zero.
inline int clear(const std::string& s) {
    const int n = static_cast<int>(s.size());
    if (n == 0) return 0;

    // flat n*n buffer, row-major on l. dp lives in one contiguous block so the
    // inner k-scan walks memory forward -- 500*500 ints is 1 MB, it stays warm.
    std::vector<int> dp(static_cast<std::size_t>(n) * n, 0);
    auto at = [&dp, n](int l, int r) -> int& {
        return dp[static_cast<std::size_t>(l) * n + r];
    };

    for (int i = 0; i < n; ++i) at(i, i) = 1;  // one letter, one operation.

    // grow the interval. every value dp reads -- dp[l+1][r], dp[l+1][k-1],
    // dp[k][r] -- spans a strictly shorter interval, so shortest-first has it
    // already computed by the time we need it.
    for (int len = 2; len <= n; ++len) {
        for (int l = 0; l + len - 1 < n; ++l) {
            const int r = l + len - 1;
            int best = at(l + 1, r) + 1;  // s[l] leaves on its own.
            for (int k = l + 1; k <= r; ++k) {
                if (s[k] != s[l]) continue;  // only equal letters can carry it.
                const int gap = (l + 1 <= k - 1) ? at(l + 1, k - 1) : 0;
                const int carried = gap + at(k, r);
                if (carried < best) best = carried;
            }
            at(l, r) = best;
        }
    }
    return at(0, n - 1);
}

}  // namespace p0006
