#pragma once

#include <string>
#include <vector>

namespace p0006 {

// the oracle. same problem, pivoted on the *last* letter instead of the first,
// written as a memoized recursion instead of a table. that's the point -- if
// the solution and the reference shared a pivot and a loop shape they'd share a
// bug, and the differential would prove nothing. this one derives itself the
// other direction so a mistake in either can't hide.
//
// dp[l][r] = cost to clear s[l..r], reasoning about how the rightmost letter
// s[r] leaves:
//   alone: dp[l][r-1] + 1
//   glued to an earlier equal s[k] (s[k] == s[r]): clear s[l..k] with s[r]
//   riding along on s[k]'s operation, then clear the gap s[k+1..r-1]:
//     dp[l][r] = dp[l][k] + dp[k+1][r-1]
//
// base: empty is zero, one letter is one.
class ReferenceSolver {
public:
    int solve(const std::string& s) {
        s_ = &s;
        const int n = static_cast<int>(s.size());
        if (n == 0) return 0;
        // -1 marks "not computed yet" -- a real cost is never negative.
        memo_.assign(static_cast<std::size_t>(n) * n, -1);
        n_ = n;
        return rec(0, n - 1);
    }

private:
    int rec(int l, int r) {
        if (l > r) return 0;
        if (l == r) return 1;
        int& m = memo_[static_cast<std::size_t>(l) * n_ + r];
        if (m >= 0) return m;

        const std::string& s = *s_;
        int best = rec(l, r - 1) + 1;  // s[r] leaves on its own.
        for (int k = l; k < r; ++k) {
            if (s[k] != s[r]) continue;
            const int cand = rec(l, k) + rec(k + 1, r - 1);
            if (cand < best) best = cand;
        }
        return m = best;
    }

    const std::string* s_ = nullptr;
    std::vector<int> memo_;
    int n_ = 0;
};

}  // namespace p0006
