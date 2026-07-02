#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0003 {

// dumb and obviously correct, for small n and q. fire every pair (i, j), paint
// the fence with everyone else, count what's lit. take the best over all pairs.
// O(q^2 * n) -- the fast solver has to agree with this on every case.
//
// no cleverness on purpose. this is the thing we trust, so it does the literal
// definition: a section is painted if some surviving painter covers it.
inline std::int64_t brute(int n, int q,
                          const std::vector<std::pair<int, int>>& ranges) {
    std::vector<char> covered(static_cast<std::size_t>(n) + 1, 0);
    std::int64_t best = 0;
    bool seen = false;

    for (int a = 0; a < q; ++a) {
        for (int b = a + 1; b < q; ++b) {
            for (int x = 1; x <= n; ++x) covered[static_cast<std::size_t>(x)] = 0;
            for (int i = 0; i < q; ++i) {
                if (i == a || i == b) continue;  // these two are fired.
                int l = ranges[static_cast<std::size_t>(i)].first;
                int r = ranges[static_cast<std::size_t>(i)].second;
                for (int x = l; x <= r; ++x) covered[static_cast<std::size_t>(x)] = 1;
            }
            std::int64_t painted = 0;
            for (int x = 1; x <= n; ++x)
                painted += covered[static_cast<std::size_t>(x)];
            if (!seen || painted > best) {
                best = painted;
                seen = true;
            }
        }
    }

    return seen ? best : 0;  // q < 2: no pair, nothing to fire.
}

}  // namespace p0003
