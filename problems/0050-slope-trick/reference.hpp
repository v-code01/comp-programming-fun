#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace p0050 {

using i64 = std::int64_t;

// dumb and obviously correct, small inputs only. the exchange argument in
// solution.hpp says an optimal b_i always lands on some input value, so restrict
// each b_i to vals = sorted distinct a. then it is a plain DP over the value grid:
//     dp_i[j] = |a_i - vals[j]| + min_{j' <= j} dp_{i-1}[j'].
// the prefix-min over j' <= j is the non-decreasing constraint spelled out --
// b_i = vals[j] may follow any earlier b_{i-1} = vals[j'] with j' <= j. the answer
// is min over j of dp_n[j]. O(n * V), V = #distinct values. no heap, no slope
// trick, no shared code with the solution -- which is the whole point of a
// differential oracle: two roads to the same number.
class ReferenceSolver {
public:
    i64 solve(const std::vector<int>& a) {
        const int n = static_cast<int>(a.size());
        if (n == 0) return 0;

        std::vector<i64> vals(a.begin(), a.end());
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        const int V = static_cast<int>(vals.size());

        // dp_0 row: cost of setting b_0 = vals[j] with no earlier constraint.
        std::vector<i64> dp(static_cast<std::size_t>(V));
        for (int j = 0; j < V; ++j) dp[j] = absdiff(a[0], vals[j]);

        std::vector<i64> ndp(static_cast<std::size_t>(V));
        const i64 INF = std::numeric_limits<i64>::max() / 4;
        for (int i = 1; i < n; ++i) {
            // running prefix-min of the previous row, then add |a_i - vals[j]|.
            i64 pref = INF;
            for (int j = 0; j < V; ++j) {
                if (dp[j] < pref) pref = dp[j];  // min over j' <= j.
                ndp[j] = pref + absdiff(a[i], vals[j]);
            }
            dp.swap(ndp);
        }
        return *std::min_element(dp.begin(), dp.end());
    }

private:
    // |x - y| in i64 -- the difference of two 1e9-scale ints overflows int32, so
    // widen before subtracting.
    static i64 absdiff(int x, i64 y) {
        i64 d = static_cast<i64>(x) - y;
        return d < 0 ? -d : d;
    }
};

}  // namespace p0050
