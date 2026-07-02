#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace p0005 {

// items weigh 1..8, cnt[i] of each. pick a subset with total <= W and make that
// total as big as you can. W is up to 1e18 and the counts up to 1e16, so the
// honest dp over [0, W] is dead on arrival -- 1e18 cells. everything below is
// about shrinking that to a few thousand and proving the shrink is exact.
//
// L = lcm(1..8) = 840. the one number the whole trick hangs on: k_i = L/i is a
// whole number for every weight, so L/i items of weight i weigh exactly L.

// bounded subset-sum reachability over [0, cap]. dp[v] == 1 iff some sub-multiset
// (a_i <= cnt_i) weighs exactly v. more than cap/i copies of weight i overshoot
// cap on their own, so cap each count there -- keeps the inner walk O(cap).
//
// the O(cap)-per-weight step is a sliding OR: along the chain v, v+i, v+2i, …,
// position k is reachable iff some earlier base dp[j]==1 sits within cnt_use
// steps behind it -- that's "add k-j copies of weight i". track the nearest base
// and you never rescan the window.
inline std::vector<char> reach(const std::array<std::int64_t, 9>& cnt,
                               std::int64_t cap) {
    std::vector<char> dp(static_cast<std::size_t>(cap) + 1, 0);
    dp[0] = 1;
    for (int i = 1; i <= 8; ++i) {
        std::int64_t cnt_use = cnt[i];
        if (cnt_use > cap / i) cnt_use = cap / i;  // extras only overshoot cap.
        if (cnt_use == 0) continue;
        std::vector<char> ndp(dp.size(), 0);
        for (std::int64_t rem = 0; rem < i && rem <= cap; ++rem) {
            std::int64_t last = std::numeric_limits<std::int64_t>::min() / 2;
            std::int64_t k = 0;
            for (std::int64_t v = rem; v <= cap; v += i, ++k) {
                if (dp[v]) last = k;  // a base reachable with zero new copies.
                if (k - last <= cnt_use) ndp[v] = 1;
            }
        }
        dp.swap(ndp);
    }
    return dp;
}

// cnt[1..8] are the counts; cnt[0] is unused. returns the max good total.
inline std::int64_t solve(std::int64_t W, const std::array<std::int64_t, 9>& cnt) {
    // total weight of everything. counts <= 1e16, weights <= 8, eight of them --
    // sum <= 36e16 = 3.6e17, well under the int64 ceiling. no overflow here.
    std::int64_t sum = 0;
    for (int i = 1; i <= 8; ++i) sum += static_cast<std::int64_t>(i) * cnt[i];

    // everything fits. take all of it -- can't do better than the whole pile.
    if (sum <= W) return sum;

    const std::int64_t L = 840;   // lcm(1..8).
    const std::int64_t TH = 16 * L;  // 13440 -- the small/large cutoff, see below.

    // --- small W: just solve it exactly. ---
    // W < 13440 means the honest dp array is tiny, so run it and stop guessing.
    // no items of weight i past floor(W/i) can appear (one alone would blow W),
    // so reach() caps them for free. the largest reachable cell <= W is the answer.
    if (W < TH) {
        std::vector<char> dp = reach(cnt, W);
        for (std::int64_t v = W; v >= 0; --v)
            if (dp[v]) return v;
        return 0;  // dp[0] is always set -- this line never actually fires.
    }

    // --- large W: residues plus full L-blocks. ---
    // sum > W here, so the answer is some achievable total just under W. every
    // achievable total T splits as T = (blocks of L) + r, r = T mod L. for a
    // fixed remainder r the *max* achievable total is what we want, then we shave
    // whole L-blocks off it until it drops under W.
    //
    // max total with remainder r = take everything, then remove the *lightest*
    // sub-multiset whose weight ≡ (sum - r) (mod L). so we need, per residue, the
    // smallest exact subset sum landing on it.
    //
    // why [0, 8L) captures every residue's lightest rep: suppose some subset of
    // weight >= 8L hits residue d. spread >= 8L across 8 weights and one weight
    // carries >= L, i.e. a_i >= k_i of it -- drop k_i of those. that's exactly -L,
    // so the residue holds and the subset stays legal. repeat while weight >= 8L.
    // you land below 8L still on residue d. so the lightest rep of any reachable
    // residue lives in [0, 8L) -- scanning that window misses nothing.
    const std::int64_t SZ = 8 * L;  // 6720.
    std::vector<char> dp = reach(cnt, SZ - 1);

    const std::int64_t INF = std::numeric_limits<std::int64_t>::max();
    std::vector<std::int64_t> min_rem(static_cast<std::size_t>(L), INF);
    for (std::int64_t v = 0; v < SZ; ++v)
        if (dp[v]) {
            std::int64_t d = v % L;
            if (v < min_rem[d]) min_rem[d] = v;  // lightest removal hitting d.
        }

    std::int64_t best = 0;  // the empty subset always fits.
    for (std::int64_t r = 0; r < L; ++r) {
        std::int64_t d = ((sum - r) % L + L) % L;  // removal must weigh ≡ d.
        if (min_rem[d] == INF) continue;           // no subset lands on r at all.
        std::int64_t mt = sum - min_rem[d];        // max achievable total ≡ r.
        if (mt < 0) continue;

        std::int64_t cand;
        if (mt <= W) {
            cand = mt;  // the max itself fits -- remove the light subset, done.
        } else {
            // shave whole L-blocks off mt until <= W. the landing is the largest
            // value ≡ r that is <= W. it exceeds W - L >= 15L, so every block we
            // remove leaves a total >= 8L -- the pigeonhole above keeps each
            // removal legal, so the landing is reachable, not just arithmetic.
            cand = W - ((W - r) % L);
        }
        if (cand > best) best = cand;
    }
    return best;
}

}  // namespace p0005
