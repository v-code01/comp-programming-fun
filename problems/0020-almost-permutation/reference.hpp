#pragma once

#include <cstdint>
#include <vector>

#include "solution.hpp"  // reuses p0020::Fact

namespace p0020 {

// the dumb oracle. tiny n only. there's nothing here to be subtly wrong about:
// tighten the intervals the same obvious way, then try every legal array a_i in
// [lo_i, hi_i] and take the smallest sum_v cnt(v)^2. -1 if any interval is empty.
//
// the sum is tracked incrementally with the same odd-number identity the flow
// leans on: dropping an item into a bucket already holding c raises cnt^2 by
// (c+1)^2 - c^2 = 2c+1. so a leaf never rescans -- each placement is O(1). this
// is the ground truth the MCMF must match, exhaustively, on thousands of cases.
//
// cost is exponential: prod of interval sizes, up to n^n. cap the differential to
// n <= 7 and it stays honest and fast.
namespace detail {

inline void rec(int i, int n, const std::vector<int>& lo,
                const std::vector<int>& hi, std::vector<int>& cnt,
                std::int64_t cur, std::int64_t& best) {
    if (i == n) {
        if (cur < best) best = cur;
        return;
    }
    for (int v = lo[static_cast<std::size_t>(i)];
         v <= hi[static_cast<std::size_t>(i)]; ++v) {
        std::int64_t add = 2 * static_cast<std::int64_t>(
                                   cnt[static_cast<std::size_t>(v)]) +
                           1;  // (c+1)^2 - c^2.
        ++cnt[static_cast<std::size_t>(v)];
        rec(i + 1, n, lo, hi, cnt, cur + add, best);
        --cnt[static_cast<std::size_t>(v)];
    }
}

}  // namespace detail

inline std::int64_t reference(int n, const std::vector<Fact>& facts) {
    std::vector<int> lo(static_cast<std::size_t>(n), 1);
    std::vector<int> hi(static_cast<std::size_t>(n), n);
    for (const Fact& f : facts) {
        for (int i = f.l; i <= f.r; ++i) {  // 1-indexed walk, map to 0-based.
            int idx = i - 1;
            if (f.t == 1) {
                if (f.v > lo[static_cast<std::size_t>(idx)]) {
                    lo[static_cast<std::size_t>(idx)] = f.v;
                }
            } else {
                if (f.v < hi[static_cast<std::size_t>(idx)]) {
                    hi[static_cast<std::size_t>(idx)] = f.v;
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (lo[static_cast<std::size_t>(i)] > hi[static_cast<std::size_t>(i)]) {
            return -1;
        }
    }

    std::vector<int> cnt(static_cast<std::size_t>(n) + 1, 0);  // cnt[1..n].
    std::int64_t best = static_cast<std::int64_t>(n) * n + 1;  // above any real cost.
    detail::rec(0, n, lo, hi, cnt, 0, best);
    return best;
}

}  // namespace p0020
