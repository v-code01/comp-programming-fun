#pragma once

#include <cstdint>
#include <queue>
#include <vector>

namespace p0050 {

using i64 = std::int64_t;

// minimum L1 cost to make an array non-decreasing. move each a[i] to any integer
// b[i], pay |a[i]-b[i]|, and land at b[1] <= b[2] <= ... <= b[n] for the least
// total. this is L1 isotonic regression. PAVA solves it. slope trick solves it.
// they are the same algorithm wearing different clothes -- the heap below IS the
// pool of the pool-adjacent-violators merge.
//
// THE CONVEX FUNCTION. sweep left to right. let
//     f_i(y) = min cost to fix the first i elements with b_i <= y.
// f_i is convex, piecewise-linear, and non-increasing -- lifting the ceiling y
// never costs more, and every piece bends the same way. the recurrence is a
// min-plus convolution:
//     f_i(y) = |a_i - y| + min_{z <= y} f_{i-1}(z).
// the inner min_{z<=y} is a prefix-minimum -- it flattens f_{i-1} to its running
// low, erasing every positive-slope piece and leaving a non-increasing envelope.
// then add |a_i - y|, a V with slope -1 then +1. min-plus with a convex |.|
// preserves convexity, so f_i stays convex. that is the invariant the whole
// thing rides on, and it is why a heap is enough: a convex PL function is nothing
// but its slope-change points.
//
// THE HEAP IS THE SLOPE-CHANGE SET. after the prefix-min flatten, f_{i-1} carries
// only non-positive slopes, so its breakpoints all sit on the LEFT arm -- the
// stretch where the slope is still negative and raising y still helps. keep
// exactly those left breakpoints in a MAX-HEAP. the top is the rightmost bend:
// the argmin ceiling, the y past which f stops falling.
//
// adding |a_i - y| drops a_i into the arm as a fresh bend. if the old top t sits
// ABOVE a_i, the V's minimum wants to pull the optimum down from t to a_i and
// that move costs t - a_i -- pay it once, drop t, seat a_i in its place:
//
//     for each a_i:
//         push a_i
//         if top() > a_i:  cost += top() - a_i;  pop();  push a_i
//
// the pop-and-repay is the entire trick and the one spot a bug hides. the "if"
// fires exactly when a_i violates monotonicity against the current optimal tail
// -- you flatten that tail down to a_i and bank the difference. nothing else.
//
// WHY b_i AT AN a-VALUE SUFFICES (the oracle leans on this). exchange argument:
// take any optimal b. if some b_i equals no a_j, slide it toward the nearest
// value it can reach without breaking b_{i-1} <= b_i <= b_{i+1}. cost |a_i - b_i|
// is linear in b_i across that free interval, so sliding to an endpoint -- an
// a-value, or a neighbor's b already pinned to an a-value -- never raises cost.
// repeat; every b_i settles on an input value. so the DP oracle that restricts
// each b_i to sorted-distinct a is exact, not an approximation.
//
// LOWER BOUND. you must read a -- Omega(n), unbeatable. slope trick sits at
// O(n log n): n pushes and at most n pops on a binary heap. the naive value-DP is
// O(n * V) over the coordinate axis; with |a_i| up to 1e9 that axis is 2e9 wide
// and the DP is dead on arrival. the log is the toll for carrying the breakpoint
// set in sorted order without paying V for the axis.
//
// int64 on the accumulator -- n up to 1e6 and each drop up to 2e9, so cost tops
// out near 2e15. values themselves fit in int32 (|a_i| <= 1e9), and int pushed
// into an i64 heap widens cleanly, so no single step overflows.
inline i64 solve(const std::vector<int>& a) {
    std::priority_queue<i64> heap;  // max-heap of left-arm breakpoints.
    i64 cost = 0;
    for (int x : a) {
        heap.push(x);
        if (heap.top() > x) {
            cost += heap.top() - x;  // flatten the tail down to x, bank the drop.
            heap.pop();
            heap.push(x);
        }
    }
    return cost;
}

}  // namespace p0050
