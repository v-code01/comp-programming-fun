#pragma once

#include <cstdint>
#include <vector>

namespace p0027 {

using i64 = std::int64_t;

// the oracle. no lagrangian, no toll, no binary search -- just the definition,
// walked day by day. at each day you either sit still, or (holding) sell, or
// (flat, budget left) buy. dp over (day, holding, trades left) is that decision
// tree with the repeats folded away. O(n*K) states, one max each. sane only for
// small n and K, which is where the differential test lives.
//
// f_i[h][j] = best profit from day i to the end, holding h positions (0 or 1),
//             with j transactions of budget remaining. a trade is spent on BUY.
//   f_i[h][j] = max( f_{i+1}[h][j],                     // skip day i
//                    h ? price_i + f_{i+1}[0][j]        // sell (must be holding)
//                      : -price_i + f_{i+1}[1][j-1] )   // buy (needs j > 0)
// answer = f_0[0][K].
inline i64 reference(int K, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    if (n == 0 || K <= 0) return 0;
    const i64 NEG = -(static_cast<i64>(1) << 62);  // an unreachable state.

    // cur = layer for day i, nxt = layer for day i+1. index [holding][budget].
    std::vector<std::vector<i64>> nxt(2, std::vector<i64>(K + 1, 0));
    std::vector<std::vector<i64>> cur(2, std::vector<i64>(K + 1, 0));

    // day n (past the end): flat with nothing owed is 0, still holding is
    // impossible to cash, so it's worth NEG -- you can't end mid-position.
    for (int j = 0; j <= K; ++j) { nxt[0][j] = 0; nxt[1][j] = NEG; }

    for (int i = n - 1; i >= 0; --i) {
        const i64 price = p[i];
        for (int j = 0; j <= K; ++j) {
            // holding: skip, or sell into the flat future.
            i64 hold = nxt[1][j];
            i64 sell = price + nxt[0][j];
            cur[1][j] = hold > sell ? hold : sell;

            // flat: skip, or (budget permitting) buy into the holding future.
            i64 flat = nxt[0][j];
            i64 buy = (j > 0) ? (nxt[1][j - 1] - price) : NEG;
            cur[0][j] = flat > buy ? flat : buy;
        }
        cur.swap(nxt);  // today becomes tomorrow for the next step down.
    }
    return nxt[0][K];  // start flat, full budget. (swap leaves day 0 in nxt.)
}

// the ground truth under the oracle: exponential, no dp, no memo. tries every
// buy/sell/skip choice literally. only for n up to ~14 -- it exists to catch a
// bug in the dp itself, so the two must never rely on the same idea.
inline i64 brute(int K, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    struct Rec {
        const std::vector<int>& p;
        int n;
        // best from day i, holding h, with j trades left. plain recursion.
        i64 go(int i, int h, int j) const {
            if (i == n) return h == 0 ? 0 : -(static_cast<i64>(1) << 62);
            i64 best = go(i + 1, h, j);  // skip.
            if (h == 1) {
                i64 s = static_cast<i64>(p[i]) + go(i + 1, 0, j);  // sell.
                if (s > best) best = s;
            } else if (j > 0) {
                i64 b = -static_cast<i64>(p[i]) + go(i + 1, 1, j - 1);  // buy.
                if (b > best) best = b;
            }
            return best;
        }
    };
    if (n == 0 || K <= 0) return 0;
    Rec rec{p, n};
    return rec.go(0, 0, K);
}

}  // namespace p0027
