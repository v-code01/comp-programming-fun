#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace p0004 {

// one charger, integer power x, you pick it once. every minute you plug it into
// exactly one laptop -- that laptop gains x and still drains its b that minute.
// a student lives if charge stays >= 0 at the start of every minute 1..k. find
// the smallest x that keeps everyone alive, or -1 if no x does.
//
// the whole thing turns on one rewrite. at the start of minute m student i holds
//   a_i - b_i*(m-1) + x * (charges it got in the first m-1 steps).
// so with s = m-1 steps of charging behind us, survival at that beginning means
//   charges_so_far >= ceil((b_i*s - a_i) / x)   whenever the top is positive.
// the right side only grows with s -- more time, more drain, more charges owed.
// so student i's demand is a staircase, and each step up is one more mandatory
// charge with a deadline. solve ceil((b_i*s - a_i)/x) >= j for the step s and the
// j-th charge lands a deadline:
//   deadline_j = floor((a_i + (j-1)*x) / b_i) + 1     (must be delivered by then)
// walk j = 1, 2, … and drop each mandatory charge (deadline <= k-1) into a bucket.
//
// now it's a scheduling question. k-1 steps, one charge each, a pile of charges
// with deadlines. earliest-deadline-first is feasible exactly when every prefix
// fits: for every step t, the charges due by t number no more than t. that's the
// classic Hall condition for unit jobs on a unit machine -- check it in one pass.
//
// x is monotone: more power shrinks every ceil and pushes every deadline later,
// so it never turns a feasible plan infeasible. binary search the smallest x.
class Solver {
public:
    // a_i, b_i in [1, 1e12]; k in [2, 2e5]. returns min x >= 0, or -1.
    std::int64_t solve(const std::vector<std::int64_t>& a,
                       const std::vector<std::int64_t>& b, std::int64_t k) {
        a_ = &a;
        b_ = &b;
        k_ = k;
        cnt_.assign(static_cast<std::size_t>(k), 0);  // deadlines land in 1..k-1.

        // if even the whole horizon's drain never sinks anyone, one charge per
        // needy student at max power is the easiest board there is -- so if this
        // is infeasible, nothing finite is. it's the search ceiling.
        //   hi = max_i b_i * (k-1) <= 1e12 * 2e5 = 2e17, well under int64's 9.2e18.
        std::int64_t maxb = 0;
        for (std::int64_t v : b) maxb = std::max(maxb, v);
        std::int64_t hi = maxb * (k - 1);
        if (!feasible(hi)) return -1;

        // smallest feasible x in [0, hi]. x = 0 is a real candidate -- feasible()
        // handles it without ever dividing by x.
        std::int64_t lo = 0;
        while (lo < hi) {
            std::int64_t mid = lo + (hi - lo) / 2;
            if (feasible(mid))
                hi = mid;
            else
                lo = mid + 1;
        }
        return lo;
    }

private:
    // can power x keep everyone alive? O(n + k): every student adds at least one
    // O(1) probe, and the total charges we ever bucket is capped at k-1 before we
    // bail -- so the staircase walks cost k across all students combined.
    bool feasible(std::int64_t x) {
        const auto& a = *a_;
        const auto& b = *b_;
        const std::int64_t n = static_cast<std::int64_t>(a.size());
        const std::int64_t steps = k_ - 1;  // charging opportunities, >= 1.

        // x = 0 charges nothing. the only survivors are the ones who never dip:
        // a_i >= b_i*(k-1). that product is <= 2e17, no overflow.
        if (x == 0) {
            for (std::int64_t i = 0; i < n; ++i)
                if (a[static_cast<std::size_t>(i)] <
                    b[static_cast<std::size_t>(i)] * steps)
                    return false;
            return true;
        }

        std::fill(cnt_.begin(), cnt_.begin() + static_cast<std::ptrdiff_t>(k_), 0);

        std::int64_t used = 0;  // mandatory charges booked so far.
        for (std::int64_t i = 0; i < n; ++i) {
            const std::int64_t bi = b[static_cast<std::size_t>(i)];
            const std::int64_t thr = bi * steps;    // <= 2e17, safe.
            std::int64_t num = a[static_cast<std::size_t>(i)];  // a_i + (j-1)*x.

            // num < thr means this charge's deadline is still inside 1..k-1, so
            // it's mandatory. once num reaches thr the student coasts to the end.
            while (num < thr) {
                std::int64_t deadline = num / bi + 1;  // in 1..k-1.
                if (++used > steps) return false;      // more charges than steps.
                cnt_[static_cast<std::size_t>(deadline)]++;
                // num < thr <= 2e17 here and x <= hi <= 2e17, so num stays under
                // 4e17 -- the add can't overflow int64.
                num += x;
            }
        }

        // earliest-deadline-first fits iff no prefix is overbooked.
        std::int64_t pref = 0;
        for (std::int64_t t = 1; t <= steps; ++t) {
            pref += cnt_[static_cast<std::size_t>(t)];
            if (pref > t) return false;
        }
        return true;
    }

    const std::vector<std::int64_t>* a_ = nullptr;
    const std::vector<std::int64_t>* b_ = nullptr;
    std::int64_t k_ = 0;
    std::vector<std::int64_t> cnt_;  // reused across the ~58 binary-search probes.
};

// thin wrapper so callers don't hold a Solver just to ask once.
inline std::int64_t min_charger_power(const std::vector<std::int64_t>& a,
                                      const std::vector<std::int64_t>& b,
                                      std::int64_t k) {
    return Solver{}.solve(a, b, k);
}

}  // namespace p0004
