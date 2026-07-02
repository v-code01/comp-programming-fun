#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace p0003 {

// integer sqrt for a known perfect square. double sqrt lands within one of the
// truth for values under 2^52, and the discriminant here is <= 2*q^2 <= 5e7 --
// nowhere near. the ±1 sweep makes it exact regardless.
inline std::int64_t isqrt(std::int64_t v) {
    if (v <= 0) return 0;
    std::int64_t d =
        static_cast<std::int64_t>(__builtin_sqrt(static_cast<double>(v)));
    while (d > 0 && d * d > v) --d;
    while ((d + 1) * (d + 1) <= v) ++d;
    return d;
}

// n sections, q painters, painter i covers [l_i, r_i]. fire exactly two of them.
// a section stays painted if any surviving painter still covers it. maximize the
// painted count.
//
// start from P = sections painted when nobody's fired. firing i and j can only
// cost you sections whose entire coverage lived inside {i, j} -- a section with
// cover >= 3 keeps a third painter who doesn't care that these two left. so:
//   loss(i, j) = only1[i] + only1[j] + only2[i][j]
//     only1[k]     = sections covered exactly once, by k alone
//     only2[i][j]  = sections covered exactly twice, by exactly {i, j}
// answer = P - min over i<j loss(i, j).
//
// getting only1 and only2 without a coverer set per section: sweep left to right
// with a difference array for the cover count, and carry two running sums over
// the active painters -- S1 = sum of ids, S2 = sum of squares. where cover is 1
// the sole id is S1 outright. where cover is 2 the pair (i, j) falls out of a
// quadratic: i + j = S1 and i^2 + j^2 = S2, so (i - j)^2 = 2*S2 - S1^2. the ids
// are < q <= 5000, everything fits int64, and the discriminant is a perfect
// square by construction -- two int64s recover the pair, no per-section storage.
//
// the min over pairs is where a naive O(q^2) scan dies for no reason. only2 is
// SPARSE -- at most n sections have cover exactly 2, so at most n pairs carry a
// nonzero only2. split the minimum:
//   forbidden pairs (only2 > 0): at most n of them, check each straight -- O(n).
//   every other pair has only2 == 0, so its loss is just only1[i] + only1[j].
//   the cheapest such pair is two small-only1 painters that aren't forbidden --
//   sort painters by only1, then for each i walk that order and take the first
//   partner that isn't i and isn't forbidden with i. i skips at most deg(i)
//   partners, and the degrees sum to 2 * (#forbidden pairs) <= 2n, so the whole
//   walk is O(n + q).
// time O(n + q log q), space O(n + q). no q*q matrix, no 25M pair scan.
inline std::int64_t solve(int n, int q,
                          const std::vector<std::pair<int, int>>& ranges) {
    const std::size_t un = static_cast<std::size_t>(n);
    std::vector<int> dcount(un + 2, 0);
    std::vector<std::int64_t> ds1(un + 2, 0);
    std::vector<std::int64_t> ds2(un + 2, 0);

    for (int i = 0; i < q; ++i) {
        int l = ranges[static_cast<std::size_t>(i)].first;
        int r = ranges[static_cast<std::size_t>(i)].second;
        std::int64_t id = i;         // 0-based -- id 0 gives S1 == 0, still fine,
        std::int64_t sq = id * id;   // the count branch decides, not S1's value.
        dcount[static_cast<std::size_t>(l)] += 1;
        dcount[static_cast<std::size_t>(r) + 1] -= 1;
        ds1[static_cast<std::size_t>(l)] += id;
        ds1[static_cast<std::size_t>(r) + 1] -= id;
        ds2[static_cast<std::size_t>(l)] += sq;
        ds2[static_cast<std::size_t>(r) + 1] -= sq;
    }

    std::vector<std::int64_t> only1(static_cast<std::size_t>(q), 0);
    // sparse only2: key i*q + j (i < j) -> count. at most n distinct keys, so
    // this is O(n) space, not O(q^2).
    std::unordered_map<std::int64_t, std::int64_t> only2;
    only2.reserve(un + 1);

    std::int64_t painted = 0;  // P: sections with cover >= 1, nobody fired.
    int count = 0;
    std::int64_t s1 = 0, s2 = 0;
    for (int x = 1; x <= n; ++x) {
        count += dcount[static_cast<std::size_t>(x)];
        s1 += ds1[static_cast<std::size_t>(x)];
        s2 += ds2[static_cast<std::size_t>(x)];
        if (count > 0) ++painted;

        if (count == 1) {
            ++only1[static_cast<std::size_t>(s1)];  // sole id is the running sum.
        } else if (count == 2) {
            std::int64_t disc = 2 * s2 - s1 * s1;  // (i - j)^2, a perfect square.
            std::int64_t d = isqrt(disc);          // = |i - j|, exact.
            int i = static_cast<int>((s1 - d) / 2);
            int j = static_cast<int>((s1 + d) / 2);
            ++only2[static_cast<std::int64_t>(i) * q + j];
        }
    }

    if (q < 2) return painted;  // can't fire two -- nothing to lose.

    std::int64_t best = -1;

    // forbidden pairs: the only ones where only2 adds to the loss.
    for (const auto& [key, cnt] : only2) {
        int i = static_cast<int>(key / q);
        int j = static_cast<int>(key % q);
        std::int64_t loss = only1[static_cast<std::size_t>(i)] +
                            only1[static_cast<std::size_t>(j)] + cnt;
        if (best < 0 || loss < best) best = loss;
    }

    // non-forbidden pairs: loss is only1[i] + only1[j]. sort by only1, then each
    // painter grabs its cheapest legal partner -- the first in order that isn't
    // itself and isn't a forbidden pair.
    std::vector<int> order(static_cast<std::size_t>(q));
    for (int i = 0; i < q; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return only1[static_cast<std::size_t>(a)] <
               only1[static_cast<std::size_t>(b)];
    });

    for (int i = 0; i < q; ++i) {
        for (int idx = 0; idx < q; ++idx) {
            int o = order[static_cast<std::size_t>(idx)];
            if (o == i) continue;
            std::int64_t key = i < o ? static_cast<std::int64_t>(i) * q + o
                                     : static_cast<std::int64_t>(o) * q + i;
            if (only2.find(key) != only2.end()) continue;  // forbidden -- skip.
            std::int64_t loss = only1[static_cast<std::size_t>(i)] +
                                only1[static_cast<std::size_t>(o)];
            if (best < 0 || loss < best) best = loss;
            break;  // first legal partner is i's cheapest -- order is ascending.
        }
    }

    if (best < 0) best = 0;
    return painted - best;
}

}  // namespace p0003
