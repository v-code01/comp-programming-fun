#pragma once

#include <cstdint>
#include <vector>

#include "solution.hpp"  // p0051::P, detail::orient

namespace p0051 {

// the face oracle. dumb, obviously correct, small n only -- the thing we trust,
// not the thing we ship.
//
// a triangle {i,j,k} is a face of the hull iff its plane SUPPORTS the point set:
// every other point lies strictly on one side of it. test it with the same exact
// predicate -- orient(i,j,k,m) over all m != i,j,k must share one nonzero sign.
// a zero would mean four coplanar points, which general position forbids; if one
// shows up the input wasn't general-position, so we discard the triple rather
// than count a degenerate face. count the triangles that pass.
//
// under general position each hull face is exactly one such triple -- three
// points on it, none other coplanar -- so the count is the number of triangular
// faces. no incremental hull, no horizon, no orientation to get wrong. just the
// definition of a supporting plane, applied C(n,3) times over n points. O(n^4).
inline int convex_hull_faces_brute(const std::vector<P>& pts) {
    using detail::orient;
    const int n = static_cast<int>(pts.size());
    if (n < 4) return 0;

    int count = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j + 1; k < n; ++k) {
                int side = 0;  // 0 = undecided, else the sign everyone shares
                bool ok = true;
                for (int m = 0; m < n; ++m) {
                    if (m == i || m == j || m == k) continue;
                    const int s = orient(pts[i], pts[j], pts[k], pts[m]);
                    if (s == 0) {  // coplanar -- not general position, not a face
                        ok = false;
                        break;
                    }
                    if (side == 0)
                        side = s;
                    else if (s != side) {  // points straddle -- not a face
                        ok = false;
                        break;
                    }
                }
                if (ok) ++count;
            }
    return count;
}

}  // namespace p0051
