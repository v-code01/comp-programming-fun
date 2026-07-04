#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "solution.hpp"  // borrows the Offer input type -- not the algorithm.

namespace p0026 {

// the independent oracle. floyd-warshall on the same -log(rate) weights, small n
// only -- it's O(n^3). different algorithm, same transform, so its verdict has to
// match bellman-ford's on every graph.
//
// build the all-pairs shortest-path table over w(u->v) = -log(rate). when it's
// done, dist[i][i] is the cheapest closed walk from i back to i. that dips below
// zero exactly when i lies on a negative cycle -- an arbitrage loop through i. so
// arbitrage exists iff some dist[i][i] < -eps.
inline bool has_arbitrage_reference(int n, const std::vector<Offer>& offers) {
    constexpr double kEps = 1e-9;
    constexpr double kInf = 1e18;

    const std::size_t sz = static_cast<std::size_t>(n) + 1;
    std::vector<std::vector<double>> dist(sz, std::vector<double>(sz, kInf));
    for (int i = 1; i <= n; ++i)
        dist[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 0.0;

    // keep the cheapest offer per (u, v). a self-loop with rate > 1 drops a
    // negative dist[u][u] right here, before the triple loop even runs.
    for (const Offer& o : offers) {
        const std::size_t u = static_cast<std::size_t>(o.u);
        const std::size_t v = static_cast<std::size_t>(o.v);
        const double weight = -std::log(o.rate);
        if (weight < dist[u][v]) dist[u][v] = weight;
    }

    for (int k = 1; k <= n; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        for (int i = 1; i <= n; ++i) {
            const std::size_t ii = static_cast<std::size_t>(i);
            if (dist[ii][kk] >= kInf * 0.5) continue;  // no i..k -- skip the row.
            for (int j = 1; j <= n; ++j) {
                const std::size_t jj = static_cast<std::size_t>(j);
                if (dist[kk][jj] >= kInf * 0.5) continue;
                const double cand = dist[ii][kk] + dist[kk][jj];
                if (cand < dist[ii][jj]) dist[ii][jj] = cand;
            }
        }
    }

    for (int i = 1; i <= n; ++i) {
        const std::size_t ii = static_cast<std::size_t>(i);
        if (dist[ii][ii] < -kEps) return true;
    }
    return false;
}

}  // namespace p0026
