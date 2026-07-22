#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "solution.hpp"  // p0054::P, detail::orient, detail::in_circle_real

namespace p0054 {

// the Delaunay-edge oracle. dumb, obviously correct, small n only -- the thing
// we trust, not the thing we ship. no triangulation is ever built.
//
// an undirected edge (i,j) is a Delaunay edge iff SOME empty circle passes
// through i and j. witness that circle by a third point k: circumcircle(i,j,k)
// is a circle through i,j, and (i,j) is Delaunay iff for SOME k it is empty of
// every other point. two directions this closes:
//   * no false negatives -- an interior Delaunay edge borders a Delaunay
//     triangle (i,j,k) with empty circumcircle; a hull edge borders exactly one
//     such triangle. either way a witness k exists.
//   * no false positives -- if circumcircle(i,j,k) is empty then (i,j,k) IS a
//     Delaunay triangle, so (i,j) is one of its edges.
// pairs x witnesses x others = O(n^4). exact via the same int128 in-circle
// predicate the solution uses -- orient the triple CCW first so its sign
// convention holds.
inline std::vector<std::pair<int, int>> delaunay_edges_brute(
    const std::vector<P>& pts) {
    using detail::in_circle_real;
    using detail::orient;
    const int n = static_cast<int>(pts.size());
    std::vector<std::pair<int, int>> edges;
    if (n < 3) return edges;

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            bool delaunay = false;
            for (int k = 0; k < n && !delaunay; ++k) {
                if (k == i || k == j) continue;
                int a = i, b = j, c = k;
                const int o = orient(pts[static_cast<std::size_t>(a)],
                                     pts[static_cast<std::size_t>(b)],
                                     pts[static_cast<std::size_t>(c)]);
                if (o == 0) continue;         // collinear -- no circle to grow
                if (o < 0) std::swap(b, c);   // make (a,b,c) CCW
                bool empty = true;
                for (int m = 0; m < n; ++m) {
                    if (m == i || m == j || m == k) continue;
                    if (in_circle_real(pts[static_cast<std::size_t>(a)],
                                       pts[static_cast<std::size_t>(b)],
                                       pts[static_cast<std::size_t>(c)],
                                       pts[static_cast<std::size_t>(m)])) {
                        empty = false;  // a site sits inside -- not this witness
                        break;
                    }
                }
                if (empty) delaunay = true;
            }
            if (delaunay) edges.emplace_back(i, j);
        }
    std::sort(edges.begin(), edges.end());  // already unique, i<j normalized
    return edges;
}

// same edge set, summed as squared lengths -- the value main prints.
inline std::int64_t delaunay_edge_sq_sum_brute(const std::vector<P>& pts) {
    std::int64_t sum = 0;
    for (const auto& e : delaunay_edges_brute(pts)) {
        const P& u = pts[static_cast<std::size_t>(e.first)];
        const P& v = pts[static_cast<std::size_t>(e.second)];
        const std::int64_t dx = u.x - v.x, dy = u.y - v.y;
        sum += dx * dx + dy * dy;
    }
    return sum;
}

}  // namespace p0054
