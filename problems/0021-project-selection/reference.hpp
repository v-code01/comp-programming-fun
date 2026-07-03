#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0021 {

// dumb and obviously correct, for tiny n only. it does not reason about cuts or
// closures -- it enumerates every one of the 2^n subsets, throws out the ones
// that break a prerequisite, and keeps the richest survivor.
//
//   a subset S (a bitmask) is CLOSED when every prerequisite it triggers is
//   satisfied: for each relation i->j, if i is in S then j must be in S too.
//   the empty set is closed and scores 0, so best starts there and the answer
//   never goes negative.
//
// exponential in n, capped at n <= 16 by the mask fitting an int -- that's the
// whole point. it's the ground truth the Dinic/min-cut solve is diffed against,
// including the nasty cases: negative profits, and prerequisite cycles that force
// a whole group all-in-or-all-out.
inline std::int64_t max_profit_brute(
    const std::vector<std::int64_t>& profits,
    const std::vector<std::pair<int, int>>& prereqs) {
    const int n = static_cast<int>(profits.size());
    std::int64_t best = 0;  // the empty closure -- always legal, always >= this.

    const unsigned limit = 1u << n;
    for (unsigned mask = 0; mask < limit; ++mask) {
        // reject the moment a triggered prerequisite is missing.
        bool closed = true;
        for (const auto& pr : prereqs) {
            const bool has_i = (mask >> pr.first) & 1u;
            const bool has_j = (mask >> pr.second) & 1u;
            if (has_i && !has_j) {
                closed = false;
                break;
            }
        }
        if (!closed) continue;

        std::int64_t sum = 0;
        for (int i = 0; i < n; ++i) {
            if ((mask >> i) & 1u) sum += profits[static_cast<std::size_t>(i)];
        }
        if (sum > best) best = sum;
    }
    return best;
}

}  // namespace p0021
