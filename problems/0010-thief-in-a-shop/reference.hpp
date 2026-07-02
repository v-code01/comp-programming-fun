#pragma once

#include <algorithm>
#include <vector>

namespace p0010 {

// the dumb one, and there's nothing in it to be subtly wrong about. the literal
// "which sums are reachable after exactly j picks":
//   reach_0 = {0};   reach_j = { s + a_i : s in reach_{j-1} }.
// run j = 1..k, read off reach_k sorted. no shift, no min-coin trick, no cap
// argument -- just the definition, walked one pick at a time over a boolean
// array. sane only for small k and small costs, which is exactly where the
// differential test feeds it.
inline std::vector<int> reference(int k, const std::vector<int>& a) {
    int maxa = 0;
    for (int x : a)
        if (x > maxa) maxa = x;

    // after j picks every sum is <= j*maxa <= k*maxa, so one array of this width
    // holds every intermediate stage without overflow.
    const int cap = k * maxa;
    std::vector<char> cur(static_cast<std::size_t>(cap) + 1, 0);
    std::vector<char> nxt(static_cast<std::size_t>(cap) + 1, 0);
    cur[0] = 1;

    for (int j = 1; j <= k; ++j) {
        std::fill(nxt.begin(), nxt.end(), 0);
        for (int s = 0; s <= cap; ++s)
            if (cur[s])
                for (int x : a)
                    nxt[s + x] = 1;  // s <= (j-1)*maxa so s + x <= j*maxa <= cap.
        cur.swap(nxt);
    }

    std::vector<int> out;
    for (int s = 0; s <= cap; ++s)
        if (cur[s]) out.push_back(s);
    return out;
}

}  // namespace p0010
