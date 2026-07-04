#pragma once

#include <cstdint>
#include <vector>

namespace p0027 {

using i64 = std::int64_t;

// optimal trading under a transaction budget. n daily prices, at most K
// round-trips (buy low, sell higher, no overlap -- one position at a time).
// maximize total profit. the naive exact dp is O(n*K); at K near n that's a
// trillion cells. the aliens trick charges a fixed toll per transaction, solves
// the toll-adjusted problem in one O(n) sweep, and binary-searches the toll to
// hit exactly K trades. O(n log maxPrice). that's the whole game.
//
// LOWER BOUND, stated plainly. you must read the prices -- Omega(n). you can't
// beat that. the trick sits at O(n log maxPrice), a log factor over the floor,
// and it buries the O(n*K) dp the moment K is large. the log is the price of not
// knowing the right toll up front.
//
// ---------------------------------------------------------------------------
// the object we optimize. g(t) = max profit using EXACTLY t transactions.
// the answer we want is max over t <= K of g(t).
//
// CONCAVITY. claim: g(t+1) - g(t) is non-increasing in t. the t-th round-trip
// you add is worth no more than the (t-1)-th. diminishing returns, provable.
//
// proof. reduce to intervals first. let d_i = p_i - p_{i-1} be the daily move.
// buying day a and selling day b collects p_b - p_a = d_{a+1} + … + d_b -- the
// sum of one contiguous block of moves. non-overlapping trades are DISJOINT
// blocks (the next buy is >= the last sell). so the problem is exactly:
//
//     pick t disjoint non-empty subarrays of d, maximize the sum of their sums.
//
// now the exchange. take an optimal (t+1)-block selection P and an optimal
// (t-1)-block selection Q. i build two selections, each with t blocks, whose
// combined weight is >= weight(P) + weight(Q). that forces
//     2*g(t) >= g(t+1) + g(t-1),
// which is concavity.
//
// lay P's blocks and Q's blocks on the line. P has two more blocks than Q. walk
// left to right and pair each Q-block with the P-block it overlaps, in order --
// a matching that leaves exactly two P-blocks unpaired (P has +2). on every
// paired stretch, cut both selections at the same coordinate and swap the tails.
// disjointness survives the swap -- the cut lands in a gap for at least one side
// at each boundary, and blocks stay blocks. weight is conserved by the swap,
// it's the same moves re-labeled. the two unpaired P-blocks are free mass; hand
// one to each side. result: two selections sized (|Q|+1) = t and (|P|-1) = t,
// total weight unchanged = weight(P) + weight(Q). each is a legal t-selection,
// so 2*g(t) >= g(t+1) + g(t-1). done.
//
// (the same fact drops out of min-cost-flow: model t trades as t units of flow,
// successive-shortest-paths augments along non-increasing profit, so the t-th
// unit earns no more than the (t-1)-th. same monotone marginal, same concavity.)
//
// concavity means g rises, plateaus, then falls -- unimodal. so max over t <= K
// of g(t) is g(K) if K sits on the rising side, and the peak value once K clears
// the peak.
// ---------------------------------------------------------------------------
//
// THE LAGRANGIAN. charge a toll lambda >= 0 per transaction. solve the
// unconstrained
//     V(lambda) = max over t of ( g(t) - lambda * t )
// in one O(n) pass -- two states, cash (flat) and hold (in a position):
//     buy  : hold <- cash - price - lambda   (open a trade, pay the toll now)
//     sell : cash <- hold + price            (close it, bank the move)
// V is the value at cash after the last day. because g is concave, g(t) -
// lambda*t is maximized at the t where the marginal g(t)-g(t-1) crosses lambda.
// raise lambda, that t drops. so lambda tunes the trade count.
//
// the conjugate identity (concave g, integer breakpoints):
//     g(K) = min over lambda >= 0 of ( V(lambda) + lambda * K ).
// weak side is free: V(lambda) >= g(K) - lambda*K for every lambda. equality
// needs a lambda whose optimal-t interval contains K.
//
// PLATEAU / TIE HANDLING -- the bug-prone spot. slopes s_t = g(t) - g(t-1) are
// integers, non-increasing. at integer lambda the set of optimal t is the
// interval { t : s_{t+1} <= lambda <= s_t }. when lambda equals some slope,
// several t tie at the same penalized value -- a plateau -- and the trade count
// jumps. so a single lambda can be optimal for a whole range of t that STRADDLES
// K. we must still read off the value at exactly K.
//
// fix: the penalized dp breaks value ties by MINIMIZING trade count. that count
// is t_lo(lambda) = |{ t : s_t > lambda }|, the low end of the plateau, and it's
// non-increasing in lambda. binary-search the SMALLEST lambda with
// t_lo(lambda) <= K. call it L. then:
//     t_lo(L)   = |{s_t > L}|      <= K        (L satisfies the search)
//     t_lo(L-1) = |{s_t > L-1}|    >= K+1  =>  |{s_t >= L}| >= K+1
//     t_hi(L)   = |{s_t >= L}|     >= K+1  >   K
// so K lands inside [t_lo(L), t_hi(L)] -- L is a valid subgradient at K even when
// the plateau jumped clean over K. equality holds:  g(K) = V(L) + L*K.
//
// we only enter this branch when K is strictly below the peak. at the peak and
// past it lambda would be 0, and a zero toll lets equal-price days spawn free
// junk trades -- so we peel the peak off separately, up front, and never let
// lambda hit 0.

// penalized sweep. returns V(lambda) and the MINIMUM trade count achieving it.
// value ties break toward fewer trades -- that minimum is t_lo(lambda), the plateau
// floor the binary search steers on.
struct Penalized {
    i64 value;  // V(lambda) = max_t ( g(t) - lambda*t )
    i64 count;  // t_lo(lambda) = fewest trades reaching that value
};

inline Penalized penalized(const std::vector<int>& p, i64 lambda) {
    // hold can't exist before the first buy -- park it below any real value.
    // -(1<<62) plus a price never overflows and never wins a max.
    const i64 NEG = -(static_cast<i64>(1) << 62);
    i64 cashV = 0, cashC = 0;    // flat: value, min trades.
    i64 holdV = NEG, holdC = 0;  // in a position: value, min trades.

    for (int price : p) {
        // both candidates read YESTERDAY's states -- so no buy-and-sell on the
        // same day (sell needs a strictly later day than its buy).
        const i64 sellV = holdV + price;         // close: bank the move.
        const i64 sellC = holdC + 1;             // one more completed trade.
        const i64 buyV = cashV - price - lambda; // open: pay price and the toll.
        const i64 buyC = cashC;                  // trade counts on close, not open.

        i64 ncV = cashV, ncC = cashC;
        if (sellV > ncV || (sellV == ncV && sellC < ncC)) { ncV = sellV; ncC = sellC; }
        i64 nhV = holdV, nhC = holdC;
        if (buyV > nhV || (buyV == nhV && buyC < nhC)) { nhV = buyV; nhC = buyC; }
        cashV = ncV; cashC = ncC;
        holdV = nhV; holdC = nhC;
    }
    // an unsold position only burns a toll -- the optimum ends flat. read cash.
    return {cashV, cashC};
}

// max profit with at most K non-overlapping buy-sell trades.
inline i64 solve(int K, const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    if (n <= 1 || K <= 0) return 0;  // no day-pair or no budget -- nothing to do.

    // the unconstrained peak U (unlimited trades) and r, the FEWEST trades to
    // reach it. U is the sum of every positive daily move. a trade can swallow a
    // zero move for free but a negative move costs profit -- so the cheapest way
    // to bank U is one trade per maximal run of non-negative moves that holds a
    // gain. r counts those runs.
    i64 U = 0;
    int r = 0;
    int maxPrice = p[0];
    bool inRun = false;  // inside a gain-carrying non-negative run.
    for (int i = 1; i < n; ++i) {
        if (p[i] > maxPrice) maxPrice = p[i];
        const int d = p[i] - p[i - 1];
        if (d < 0) {
            inRun = false;                 // a loss move breaks the run.
        } else if (d > 0) {
            U += d;                         // bank the gain.
            if (!inRun) { ++r; inRun = true; }  // a new run opens one trade.
        }
        // d == 0 rides the current run -- free to hold across, opens nothing.
    }

    // budget clears the cheapest optimum -- take the whole peak, no toll needed.
    // this also swallows the K >= n/2 "unlimited" case and the all-down / flat
    // cases (r == 0, U == 0).
    if (K >= r) return U;

    // K < r: strictly on the rising side, g(K) < U. tune the toll.
    // lambda* lies in [1, best single-trade profit] <= maxPrice. it's never 0
    // here -- that's exactly what the K >= r peel guaranteed.
    i64 lo = 1, hi = maxPrice;
    while (lo < hi) {
        const i64 mid = lo + (hi - lo) / 2;
        if (penalized(p, mid).count <= K) hi = mid;  // toll high enough -- pull down.
        else lo = mid + 1;                            // too many trades -- charge more.
    }
    const Penalized best = penalized(p, lo);
    // L = lo is a valid subgradient at K (plateau argument above), so this is
    // exactly g(K).
    return best.value + lo * static_cast<i64>(K);
}

}  // namespace p0027
