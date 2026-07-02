#pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <vector>

namespace p0010 {

// take EXACTLY k products, repetition allowed, prices a_i. which grand totals
// can you land on? the whole thing turns on a shift. let m = min a_i. every pick
// costs at least m, so any k-pick total is k*m plus a non-negative "extra" s --
// the sum of (a_i - m) over the k items. set b_i = a_i - m. the cheapest kind
// has b = 0. that zero is the load-bearing part.
//
// claim: total T is reachable with EXACTLY k picks
//        iff  T = k*m + s  and  dp[s] <= k,
// where dp[s] = fewest items (unbounded supply) whose b-values sum to exactly s.
//
// proof, both directions --
//   (<=) dp[s] = t <= k means some t items hit b-sum s. a zero-b item exists,
//        the cheapest kind, so toss in k - t copies of it. now it's exactly k
//        picks, the b-sum is still s (zeros move nothing), and the a-cost is
//        s + k*m. this is the entire "pad with zeros" step: t <= k is precisely
//        what lets you inflate a t-item multiset to k items without touching the
//        sum. no zero-cost kind, no padding -- but m itself is that kind.
//   (=>) k picks with a-sum T give b-sum s = T - k*m >= 0. those k items reach s,
//        so the minimum to reach s is at most k. dp[s] <= k. done.
//
// so the answer set is exactly { k*m + s : dp[s] <= k }, and dp is coin-change
// min-coins over the DISTINCT POSITIVE b-values -- a zero coin never lowers a
// minimum, the padding already owns it. s never exceeds k*maxb: k items cap the
// b-sum right there, so the array stops at k*maxb and misses nothing above it.
//
// why this, not fft -- V = k*maxb <= ~1e6, D = distinct positive costs <= ~1e3.
// the dp is O(V*D) ~ 1e9 tight integer ops, one load and one min per coin, over
// a back-window of <= maxb that stays hot in cache. exact arithmetic, no
// rounding. an fft k-th power of the support polynomial is O(V log V log k) on
// paper -- but exactly-k path counts blow past double range fast, so you end up
// thresholding float noise to recover the support, and that's a real correctness
// hazard. at these bounds the integer dp wins on both the constant factor and the
// blast radius. it's also the accepted solution.
//
// the honest floor: the output alone can carry ~k*maxb ~ 1e6 values, so any
// solver is Omega(V) just to emit, plus Omega(n) to read the costs. this does
// O(n + V*D) -- within the coin count D of that floor, and you can't beat
// Omega(V) when the answer itself is that wide.
//
// returns every achievable total, strictly increasing.
inline std::vector<int> solve(int k, const std::vector<int>& a) {
    // one pass for the min and the max. n >= 1, so a[0] is a valid seed.
    int m = a[0];
    int maxa = a[0];
    for (int x : a) {
        if (x < m) m = x;
        if (x > maxa) maxa = x;
    }

    // the coins: distinct positive shifts b_i = a_i - m. sorted, so the inner
    // loop can bail the moment a coin outgrows the target.
    std::vector<int> coins;
    coins.reserve(a.size());
    for (int x : a)
        if (x - m > 0) coins.push_back(x - m);
    std::sort(coins.begin(), coins.end());
    coins.erase(std::unique(coins.begin(), coins.end()), coins.end());

    const int maxb = maxa - m;          // 0 when every price is equal.
    const int V = k * maxb;             // <= 1000*999 < 1e6 -- fits int easily.
    const int INF = INT_MAX / 2;        // sentinel above any real coin count.

    // dp[s] = min coins summing to exactly s, capped at INF for unreachable s.
    std::vector<int> dp(static_cast<std::size_t>(V) + 1, INF);
    dp[0] = 0;                          // zero coins make zero.
    for (int s = 1; s <= V; ++s) {
        int best = INF;
        for (int c : coins) {
            if (c > s) break;           // sorted -- everything past here is bigger.
            int prev = dp[s - c];       // s - c >= 0, already computed.
            if (prev < best) best = prev;
        }
        dp[s] = (best >= INF) ? INF : best + 1;  // one more coin on the best base.
    }

    // emit k*m + s for every extra-sum reachable within k items. ascending in s
    // means ascending in total, so no final sort.
    std::vector<int> out;
    for (int s = 0; s <= V; ++s)
        if (dp[s] <= k) out.push_back(k * m + s);
    return out;
}

}  // namespace p0010
