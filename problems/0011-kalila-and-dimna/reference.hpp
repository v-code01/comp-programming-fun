#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0011 {

// the dumb, obviously-correct oracle. small inputs only.
//
// it is the recurrence, transcribed with nothing clever bolted on:
//   dp[0] = 0                                    -- first tree (height 1) is free.
//   dp[i] = min over j < i of dp[j] + a_i*b_j.
// answer = dp[n-1]. no hull, no envelope, no monotonicity assumed -- every i scans
// every j < i outright. O(n^2), which is exactly why it stays trustworthy: there
// is no line geometry to get the sign of wrong. the fast solution has to match it
// on thousands of random valid boards.
//
// same overflow reasoning as the real one: a_i*b_j <= 1e18 and dp[i] <= a_i*b_1
// <= 1e18, so nothing here wraps int64 -- but the tests keep values modest anyway,
// so a silent wrap would surface as a diff rather than hide.
inline std::int64_t reference_min_cost(const std::vector<std::int64_t>& a,
                                       const std::vector<std::int64_t>& b) {
    const std::size_t n = a.size();
    std::vector<std::int64_t> dp(n, 0);
    for (std::size_t i = 1; i < n; ++i) {
        std::int64_t best = dp[0] + a[i] * b[0];  // j = 0 always available.
        for (std::size_t j = 1; j < i; ++j) {
            const std::int64_t cand = dp[j] + a[i] * b[j];
            if (cand < best) best = cand;
        }
        dp[i] = best;
    }
    return dp[n - 1];
}

}  // namespace p0011
