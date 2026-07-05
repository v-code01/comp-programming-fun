#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

namespace p0032 {

// the dumb, obviously-correct oracle -- subset DP over vertex masks, exact for
// tiny n. it shares NOT ONE line with the blossom: no duals, no trees, no
// shrink/expand. it just asks, for every subset of vertices, "what's the best
// matching that uses only these?" and answers by looking at the lowest vertex.
//
//   dp[mask] = max weight of a matching drawn from the vertices in mask.
//   take i = the lowest set vertex of mask. in any matching i is either
//     unmatched  -> dp[mask without i], or
//     matched to some j in mask via edge (i, j) -> w(i,j) + dp[mask without i,j].
//   the best over "leave it" and every legal partner j is dp[mask]. dp[0] = 0,
//   so leaving everyone unmatched is always on the table and the answer is >= 0.
//
// O(2^n * n), capped at n <= 16 by the mask fitting an unsigned. this IS the
// ground truth the blossom is diffed against -- especially odd cycles, where the
// blossom actually forms and the bipartite shortcuts would lie.
inline std::int64_t max_weight_matching_ref(
    int n, const std::vector<std::tuple<int, int, std::int64_t>>& edges) {
    if (n <= 0) return 0;
    const std::size_t N = static_cast<std::size_t>(n);

    // weight matrix, deduped by max. absent / non-positive stays 0 and never
    // gets picked (a 0 edge can't beat leaving the vertex free).
    std::vector<std::vector<std::int64_t>> w(N, std::vector<std::int64_t>(N, 0));
    for (const auto& [a, b, ww] : edges) {
        if (a == b || ww <= 0) continue;
        if (ww > w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]) {
            w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = ww;
            w[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] = ww;
        }
    }

    const unsigned full = 1u << n;
    std::vector<std::int64_t> dp(full, 0);
    for (unsigned mask = 1; mask < full; ++mask) {
        int i = __builtin_ctz(mask);  // lowest set vertex -- decide its fate.
        std::int64_t best = dp[mask & ~(1u << i)];  // leave i unmatched.
        for (int j = i + 1; j < n; ++j) {
            if (!((mask >> j) & 1u)) continue;
            const std::int64_t e =
                w[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            if (e <= 0) continue;  // no edge -- can't pair here.
            const std::int64_t cand =
                e + dp[mask & ~(1u << i) & ~(1u << j)];
            if (cand > best) best = cand;
        }
        dp[mask] = best;
    }
    return dp[full - 1];
}

}  // namespace p0032
