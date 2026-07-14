#pragma once

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

#include "solution.hpp"  // Edge, kInf -- the types, none of the reasoning.

namespace p0041 {

// dumb and obviously correct, tiny n only. no subset DP over terminals anywhere
// -- if the oracle shared the idea it's checking, it would share the bug.
//
// the whole argument in one line: an optimal Steiner tree is a spanning tree of
// its OWN vertex set. so enumerate that vertex set directly.
//
//   for every subset S of the n vertices that contains all terminals:
//       if the induced subgraph on S is connected,
//       take the weight of its minimum spanning tree.
//   answer = the min over all such S.
//
// why it's right in both directions:
//   >=  every candidate here is a connected subgraph containing the terminals,
//       so the true optimum is no worse than the best of them.
//   <=  let T be an optimal Steiner tree and S its vertex set. S contains the
//       terminals, the induced subgraph on S is connected (T alone connects it),
//       and T is a spanning tree of S -- so MST(S) <= w(T). the enumeration sees
//       S, so it reports at most w(T).
//
// Kruskal does double duty: it returns the MST weight and tells us whether the
// induced subgraph was connected at all (one component left, or nothing doing).
//
// O(2^n * m alpha) with the edges sorted once up front. n <= ~14 in tests --
// 2^n is the point, and the point is that it can't be cheating.
inline std::int64_t steiner_tree_brute(int n, const std::vector<Edge>& edges,
                                       const std::vector<int>& terminals) {
    if (terminals.size() <= 1) return 0;

    unsigned need = 0;
    for (int t : terminals) need |= 1u << t;

    // sort once, outside the 2^n loop. Kruskal only ever wants the order.
    std::vector<Edge> sorted = edges;
    std::sort(sorted.begin(), sorted.end(),
              [](const Edge& a, const Edge& b) { return a.w < b.w; });

    std::vector<int> parent(static_cast<std::size_t>(n));
    std::int64_t best = kInf;

    const unsigned limit = 1u << n;
    for (unsigned S = 0; S < limit; ++S) {
        if ((S & need) != need) continue;  // must carry every terminal.

        std::iota(parent.begin(), parent.end(), 0);
        int comps = __builtin_popcount(S);  // vertices in S, each its own island.
        std::int64_t total = 0;

        // plain union-find with path compression, written inline and flat.
        auto find = [&parent](int x) {
            while (parent[static_cast<std::size_t>(x)] != x) {
                parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(
                    parent[static_cast<std::size_t>(x)])];
                x = parent[static_cast<std::size_t>(x)];
            }
            return x;
        };

        for (const Edge& e : sorted) {
            if (comps == 1) break;  // spanning already. the rest is heavier.
            if (!((S >> e.u) & 1u) || !((S >> e.v) & 1u)) continue;  // not induced.
            const int ru = find(e.u);
            const int rv = find(e.v);
            if (ru == rv) continue;  // would close a cycle.
            parent[static_cast<std::size_t>(ru)] = rv;
            total += e.w;
            --comps;
        }

        if (comps != 1) continue;      // induced subgraph is disconnected. skip.
        if (total < best) best = total;
    }
    return best;
}

}  // namespace p0041
