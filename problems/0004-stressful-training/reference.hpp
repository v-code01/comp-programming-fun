#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace p0004 {

// the dumb, obviously-correct oracle. small inputs only.
//
// feasibility by direct simulation: carry every laptop's live charge, and at each
// step hand the charger to whoever dies soonest -- the student with the fewest
// whole minutes of runway left, floor(charge / drain). that greedy is exactly
// right: if an optimal plan ever charges someone else while the soonest-to-die
// still needs help, swap the two charges. the swap can't hurt the one who had
// slack (they had time to spare) and it saves the one who didn't. so charging the
// most urgent student first never loses a feasible board. extra charge on a
// student who's already safe only raises their charge -- it never sinks anyone.
// so "always charge the most urgent" decides feasibility.
//
// min-x finder: x = 0, 1, 2, … up to the ceiling hi = max_i b_i*(k-1). first x
// that survives is the answer; none by hi means -1. no monotonicity assumed and
// no binary search -- it checks each x outright, which is the point of an oracle.
// kept honest by keeping the test inputs tiny, so hi stays small.
class ReferenceSolver {
public:
    std::int64_t solve(const std::vector<std::int64_t>& a,
                       const std::vector<std::int64_t>& b, std::int64_t k) {
        std::int64_t maxb = 0;
        for (std::int64_t v : b) maxb = std::max(maxb, v);
        std::int64_t hi = maxb * (k - 1);
        for (std::int64_t x = 0; x <= hi; ++x)
            if (feasible(a, b, k, x)) return x;
        return -1;
    }

private:
    static bool feasible(const std::vector<std::int64_t>& a,
                         const std::vector<std::int64_t>& b, std::int64_t k,
                         std::int64_t x) {
        std::vector<std::int64_t> e = a;  // charge at the start of the minute.
        const std::size_t n = a.size();

        // k-1 charging steps, one per minute. after the step, drain everyone and
        // check the next beginning is still >= 0.
        for (std::int64_t t = 1; t <= k - 1; ++t) {
            // most urgent = smallest whole-minute runway left.
            std::size_t c = 0;
            std::int64_t best = e[0] / b[0];
            for (std::size_t i = 1; i < n; ++i) {
                std::int64_t runway = e[i] / b[i];
                if (runway < best) {
                    best = runway;
                    c = i;
                }
            }

            for (std::size_t i = 0; i < n; ++i) e[i] -= b[i];  // this minute's drain.
            e[c] += x;                                          // the one we charged.

            for (std::size_t i = 0; i < n; ++i)
                if (e[i] < 0) return false;  // someone hit the next minute negative.
        }
        return true;
    }
};

}  // namespace p0004
