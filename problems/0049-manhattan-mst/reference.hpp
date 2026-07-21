#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "solution.hpp"  // for p0049::Point

namespace p0049 {

// the oracle. no octants, no fenwick, no cleverness -- the complete graph and
// dense prim, straight out of the textbook. for every point not yet in the tree
// keep its cheapest edge to the tree; each round pull in the cheapest such point
// and relax the rest. O(n^2). it looks at ALL C(n,2) manhattan distances, so if
// the fast solution ever skips a candidate it should skip, this catches it.
//
// small n only -- this is the thing we trust, not the thing we ship.
class ReferenceSolver {
public:
    std::int64_t solve(const std::vector<Point>& pts) {
        const int n = static_cast<int>(pts.size());
        if (n <= 1) return 0;

        constexpr std::int64_t kInf = std::numeric_limits<std::int64_t>::max();
        std::vector<std::int64_t> dist(n, kInf);  // cheapest edge into the tree
        std::vector<char> in_tree(n, 0);
        dist[0] = 0;

        std::int64_t total = 0;
        for (int iter = 0; iter < n; ++iter) {
            // pick the closest point still outside the tree.
            int u = -1;
            std::int64_t best = kInf;
            for (int v = 0; v < n; ++v) {
                if (!in_tree[v] && dist[v] < best) {
                    best = dist[v];
                    u = v;
                }
            }
            // n is small and the graph is complete, so u is always found.
            in_tree[u] = 1;
            total += dist[u];

            // relax every other point against the newcomer.
            for (int v = 0; v < n; ++v) {
                if (in_tree[v]) continue;
                std::int64_t d = manhattan(pts[u], pts[v]);
                if (d < dist[v]) dist[v] = d;
            }
        }
        return total;
    }

private:
    static std::int64_t manhattan(const Point& a, const Point& b) {
        std::int64_t dx = a.x - b.x, dy = a.y - b.y;
        return (dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy);
    }
};

}  // namespace p0049
