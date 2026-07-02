#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0011 {

// n trees. tree i has height a_i and a recharge cost b_i. every chop drops one
// tree by 1, and between chops the saw has to recharge -- the recharge costs
// b_j where j is the tallest-indexed tree already flattened to 0. cut them all.
// minimize the total recharge bill.
//
// the guarantees do the heavy lifting: a strictly increasing, b strictly
// decreasing, a_1 = 1, b_n = 0. a_1 = 1 means tree 1 falls in one chop off the
// saw's initial charge -- so felling tree 1 is free. b_n = 0 means once the
// tallest tree is down every later recharge costs nothing -- so the job ends the
// instant tree n hits 0. what's left to price is the recharges you buy on the way
// to felling each tree, and to chop tree i a_i times you want some already-felled
// tree j whose per-recharge cost b_j you pay once per chop.
//
// dp[i] = min total recharge cost to fully fell tree i.
//   dp[1] = 0                                  -- free off the initial charge.
//   dp[i] = min over j < i of dp[j] + a_i*b_j  -- fell some j first (cost dp[j]),
//                                                 then pay b_j on each of the a_i
//                                                 chops tree i takes.
// answer = dp[n]. j < i is enough: b is strictly decreasing so any j you'd want
// to lean on -- a small per-chop cost -- has a smaller index than i anyway.
//
// the recurrence is a line-min query wearing a disguise. carry each solved j as a
// line
//   f_j(x) = b_j * x + dp[j]        -- slope b_j, intercept dp[j].
// then dp[i] = min_j f_j(a_i): evaluate every line at x = a_i, keep the lowest.
// that's the CONVEX HULL TRICK -- store only the lines on the lower envelope and
// each query is one walk to the winning line.
//
// two monotonicities collapse it to a deque with no binary search anywhere:
//   slopes b_j are STRICTLY DECREASING in j, and we insert j in index order, so
//     every new line's slope undercuts every line present -- it always appends to
//     the back of the hull. no sorted insert.
//   queries a_i are STRICTLY INCREASING in i. on a lower envelope of lines the
//     winning slope only falls as x climbs, so the answer pointer only walks
//     forward, never back. amortized O(1) per insert and per query -> O(n) total.
//
// that O(n) is optimal to a constant: you can't answer without reading all 2n
// numbers, so Omega(n) is a hard floor and the hull sits right on it. a Li Chao
// tree would also solve this at O(n log(max a)) -- strictly worse here because the
// monotone structure is handed to us for free.
//
// overflow: a_i, b_j <= 1e9 so a_i*b_j <= 1e18, and dp[i] = min of single such
// terms so dp[i] <= a_i*b_1 <= 1e18 -- a line's value m*x + c stays under 2e18,
// well inside int64's 9.2e18. the one place int64 is NOT enough is the hull's
// cross-product test, called out at bad() below.
class HullMin {
public:
    void reserve(std::size_t n) {
        m_.reserve(n);
        c_.reserve(n);
    }

    // append line y = m*x + c. slopes arrive strictly decreasing (b is strictly
    // decreasing), so the new line can only retire lines at the BACK -- pop each
    // one it hides, then push.
    void add(std::int64_t m, std::int64_t c) {
        while (m_.size() >= 2 && bad(m_.size() - 2, m_.size() - 1, m, c)) {
            m_.pop_back();
            c_.pop_back();
        }
        m_.push_back(m);
        c_.push_back(c);
    }

    // min over all stored lines at x. queries arrive non-decreasing, so the
    // pointer only advances -- amortized O(1) across the whole run. the clamp
    // catches the case where a pop retired the line the pointer was resting on;
    // when that happens the newly pushed line already dominates at this x and
    // every later (larger) x, so resuming from the back and walking forward lands
    // on the true winner.
    std::int64_t query(std::int64_t x) {
        if (ptr_ >= m_.size()) ptr_ = m_.size() - 1;
        while (ptr_ + 1 < m_.size() && eval(ptr_ + 1, x) <= eval(ptr_, x)) ++ptr_;
        return eval(ptr_, x);
    }

private:
    // m <= 1e9 and x = a_i <= 1e9 so m*x <= 1e18; c = dp <= 1e18. the sum stays
    // under 2e18 < 9.2e18 -- int64 holds every line value with room to spare.
    std::int64_t eval(std::size_t i, std::int64_t x) const {
        return m_[i] * x + c_[i];
    }

    // is line i2 (the current back) dead once the line (m3, c3) joins? the three
    // lines i1, i2, (m3,c3) have strictly decreasing slopes. on a lower envelope
    // line i2 only surfaces if its takeover point sits strictly left of the next
    // line's -- intersect(i1,i2) < intersect(i2,new). when that fails i2 is never
    // the minimum for any x, so drop it.
    //   intersect(p,q) = (c_q - c_p) / (m_p - m_q),  m_p > m_q so the denom > 0.
    // clearing both positive denominators (which flips nothing), i2 is redundant iff
    //   (c2 - c1)*(m2 - m3) >= (c3 - c2)*(m1 - m2).
    // THE OVERFLOW TRAP: those cross terms blow past int64. (c2 - c1) reaches ~1e18
    // and (m2 - m3) reaches ~1e9, so the product reaches ~1e27 -- 100 million times
    // int64's ceiling. compute it in __int128 (ceiling ~1.7e38) where it cannot
    // wrap; the compared magnitudes are exact and the >= is honest.
    bool bad(std::size_t i1, std::size_t i2, std::int64_t m3,
             std::int64_t c3) const {
        const std::int64_t m1 = m_[i1], c1 = c_[i1];
        const std::int64_t m2 = m_[i2], c2 = c_[i2];
        const __int128 lhs = static_cast<__int128>(c2 - c1) * (m2 - m3);
        const __int128 rhs = static_cast<__int128>(c3 - c2) * (m1 - m2);
        return lhs >= rhs;
    }

    std::vector<std::int64_t> m_;  // slopes, strictly decreasing on the hull.
    std::vector<std::int64_t> c_;  // matching intercepts.
    std::size_t ptr_ = 0;          // last winner; monotone queries never rewind it.
};

// a_i strictly increasing with a_1 = 1, b_i strictly decreasing with b_n = 0.
// returns the minimum total recharge cost to fell every tree. O(n) time, O(n)
// space. n == 1 returns 0 -- the lone tree is height 1, felled off the initial
// charge for free.
inline std::int64_t min_cost(const std::vector<std::int64_t>& a,
                             const std::vector<std::int64_t>& b) {
    const std::size_t n = a.size();
    HullMin hull;
    hull.reserve(n);

    std::int64_t dp = 0;      // dp for the first tree (0-indexed): free.
    hull.add(b[0], dp);       // its line, ready for later queries.
    for (std::size_t i = 1; i < n; ++i) {
        dp = hull.query(a[i]);  // min over j < i of dp[j] + a_i*b_j.
        hull.add(b[i], dp);     // this tree becomes a candidate for the rest.
    }
    return dp;                  // dp for the last tree = dp[n].
}

}  // namespace p0011
