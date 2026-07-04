#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace p0028 {

// the oracle. no suffix array, no lcp, no sorting -- the textbook O(|s|*|t|) dp.
// this is what the fast solution has to agree with on every small pair, and it
// shares not one line of machinery with it. if both were suffix structures a bug
// could live in both; this one can't hide there.
//
// dp[i][j] = length of the common substring that *ends* at s[i-1] and t[j-1].
// the two characters either match and extend the diagonal run, or they don't and
// the run resets to zero:
//   s[i-1] == t[j-1]  ->  dp[i][j] = dp[i-1][j-1] + 1
//   else              ->  dp[i][j] = 0
// the answer is the largest cell -- the longest run that ended anywhere.
//
// only the previous row is ever read, so two rolling rows hold it in O(|t|)
// space. small inputs only; this exists to be obviously right, not fast.
inline int reference_lcs(const std::string& s, const std::string& t) {
    const int n = static_cast<int>(s.size());
    const int m = static_cast<int>(t.size());
    if (n == 0 || m == 0) return 0;

    std::vector<int> prev(m + 1, 0), cur(m + 1, 0);
    int best = 0;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (s[i - 1] == t[j - 1]) {
                cur[j] = prev[j - 1] + 1;
                if (cur[j] > best) best = cur[j];
            } else {
                cur[j] = 0;
            }
        }
        std::swap(prev, cur);
    }
    return best;
}

}  // namespace p0028
