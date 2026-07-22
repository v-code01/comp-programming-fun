#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0051 {

// 3D convex hull -- count the triangular faces.
//
// n points in space, integer coords, general position: no 3 collinear, no 4
// coplanar. so the hull is simplicial -- every face is a triangle -- and the
// job is just to report how many.
//
// -- the count is pinned by Euler --
// V - E + F = 2 for any convex polytope. simplicial means every face has 3
// edges and every edge borders 2 faces, so 3F = 2E. drop E = 3F/2 in:
//   V - 3F/2 + F = 2   =>   F = 2V - 4.
// V here is the number of points that actually land on the hull -- interior
// points don't count toward it. so the answer is 2*(hull vertices) - 4, and it
// equals faces.size() once the surface is built. we build the surface and
// count. no shortcut -- we don't know V up front.
//
// -- the predicate everything rides on --
// one sign decides which side of an oriented plane a point sits.
//   orient(a,b,c,d) = sign of the signed volume of tetra (a,b,c,d)
//                   = sign of det[b-a; c-a; d-a]
//                   = sign of (b-a) x (c-a) . (d-a).
// > 0 : d is on the side the normal (b-a)x(c-a) points to -- "in front".
// < 0 : behind.   = 0 : coplanar -- general position rules it out.
//
// coords reach 1e6, so a difference reaches 2e6. the cross product squares that
// (~8e12), the dot multiplies once more (~2e6) and sums three terms -- the
// triple product reaches ~5e19. int64 caps at 9.2e18. it WOULD overflow. so the
// product runs in __int128. no float anywhere -- the sign has to be exact or the
// hull tears along a face it read the wrong way.
//
// -- incremental construction, O(n^2) --
// seed a tetra from the first 4 points -- non-coplanar, guaranteed. orient each
// face outward: its opposite vertex sits behind it. outward normals on a convex
// body induce a consistent orientation -- every edge gets walked once in each
// direction -- and that consistency is the whole trick below.
//
// then add the rest one at a time. for a new point p:
//   a face is VISIBLE if p is in front of it, orient(face, p) > 0 -- p sees it
//   from outside, so the face can't survive: p is beyond it.
//   the visible faces form one connected patch. delete them. the boundary of
//   that patch -- the HORIZON -- is a cycle of edges shared with kept faces. fan
//   p to every horizon edge with a new triangle. that's the updated hull.
//   if p sees nothing, p is inside. drop it -- it never reaches the surface.
//
// -- the horizon, without adjacency pointers --
// each face stores its vertices in a fixed cyclic order (a,b,c), so it owns three
// DIRECTED edges a->b, b->c, c->a. consistent orientation means an interior edge
// {u,v} shows up as u->v in one face and v->u in the other. so:
//   stamp every directed edge of the visible faces.
//   a directed edge u->v is on the horizon iff its reverse v->u is NOT stamped
//   -- the face across it is being kept.
// each such u->v becomes a new face (u, v, p). why that order: the kept neighbor
// walks the shared edge v->u, so the new face must walk it u->v to stay
// consistent -- neighbors traverse a shared edge in opposite directions. the
// same argument fans the new triangles consistently around p. orientation is
// inherited from the surface, never recomputed. no flips, no centroid, no float.
//
// -- cost, and the floor --
// reading the points is Omega(n). each insertion scans every current face for
// visibility: O(F) = O(n) faces, since F = 2V-4. n insertions -> O(n^2). a
// randomized-incremental or divide-and-conquer hull hits O(n log n); we take
// O(n^2) on purpose -- at n <= 2000 that's ~8e6 predicate calls, and the flat
// loop is easier to prove correct than a conflict graph. exactness over
// asymptotics.

struct P {
    std::int64_t x, y, z;
};

namespace detail {

// sign of the signed volume of tetra (a,b,c,d): +1 / 0 / -1. the __int128 is
// load-bearing -- the triple product reaches ~5e19, past int64.
inline int orient(const P& a, const P& b, const P& c, const P& d) {
    const __int128 ux = (__int128)b.x - a.x, uy = (__int128)b.y - a.y,
                   uz = (__int128)b.z - a.z;
    const __int128 vx = (__int128)c.x - a.x, vy = (__int128)c.y - a.y,
                   vz = (__int128)c.z - a.z;
    const __int128 wx = (__int128)d.x - a.x, wy = (__int128)d.y - a.y,
                   wz = (__int128)d.z - a.z;
    const __int128 vol = ux * (vy * wz - vz * wy)   //
                         - uy * (vx * wz - vz * wx)  //
                         + uz * (vx * wy - vy * wx);
    return (vol > 0) - (vol < 0);
}

struct Face {
    int a, b, c;
};

}  // namespace detail

// build the hull of pts and return its triangular face count. wants n >= 4
// points in general position -- the problem guarantees it. O(n^2).
inline int convex_hull_faces(const std::vector<P>& pts) {
    using detail::Face;
    using detail::orient;
    const int n = static_cast<int>(pts.size());
    // fewer than 4 points bound no volume -- no faces. the statement gives
    // n >= 4; the guard keeps the oracle differential honest on tiny inputs.
    if (n < 4) return 0;

    // seed the tetra from 0,1,2,3. general position => non-coplanar => a real
    // tetra. orient each face so its opposite vertex sits behind it (outward).
    std::vector<Face> faces;
    faces.reserve(static_cast<std::size_t>(2) * n);
    auto add_tetra_face = [&](int a, int b, int c, int opp) {
        // opp behind the face means orient(a,b,c,opp) < 0. if it's in front,
        // swap two vertices -- that flips the normal to point outward.
        if (orient(pts[a], pts[b], pts[c], pts[opp]) > 0) std::swap(b, c);
        faces.push_back({a, b, c});
    };
    add_tetra_face(0, 1, 2, 3);
    add_tetra_face(0, 1, 3, 2);
    add_tetra_face(0, 2, 3, 1);
    add_tetra_face(1, 2, 3, 0);

    // horizon bookkeeping. stamp[u*n + v] == epoch  <=>  u->v is a directed edge
    // of a visible face this round. a flat array plus an epoch counter beats
    // clearing n^2 slots every insertion -- bump the epoch and the old stamps go
    // stale for free. n <= 2000 => 4e6 ints, ~16MB, allocated once.
    std::vector<int> stamp(static_cast<std::size_t>(n) * n, 0);
    int epoch = 0;

    std::vector<char> vis;      // per-face visibility, reused each round
    std::vector<Face> horizon;  // the fan of new faces this round
    std::vector<Face> next;     // rebuilt face list this round

    for (int p = 4; p < n; ++p) {
        ++epoch;
        const int fcount = static_cast<int>(faces.size());
        vis.assign(static_cast<std::size_t>(fcount), 0);

        bool any = false;
        for (int i = 0; i < fcount; ++i) {
            const Face& f = faces[static_cast<std::size_t>(i)];
            if (orient(pts[f.a], pts[f.b], pts[f.c], pts[p]) > 0) {
                vis[static_cast<std::size_t>(i)] = 1;
                any = true;
                // stamp the face's three directed edges.
                stamp[static_cast<std::size_t>(f.a) * n + f.b] = epoch;
                stamp[static_cast<std::size_t>(f.b) * n + f.c] = epoch;
                stamp[static_cast<std::size_t>(f.c) * n + f.a] = epoch;
            }
        }
        if (!any) continue;  // p sits inside the hull -- it never surfaces.

        // horizon = visible directed edges whose reverse is not visible. the
        // reverse being stamped means the neighbor is also going away, so that
        // edge is interior to the deleted patch -- skip it.
        horizon.clear();
        auto emit = [&](int u, int v) {
            if (stamp[static_cast<std::size_t>(v) * n + u] != epoch)
                horizon.push_back({u, v, p});
        };
        for (int i = 0; i < fcount; ++i) {
            if (!vis[static_cast<std::size_t>(i)]) continue;
            const Face& f = faces[static_cast<std::size_t>(i)];
            emit(f.a, f.b);
            emit(f.b, f.c);
            emit(f.c, f.a);
        }

        // drop visible faces, keep the rest, append the fan around p.
        next.clear();
        next.reserve(static_cast<std::size_t>(fcount) + horizon.size());
        for (int i = 0; i < fcount; ++i)
            if (!vis[static_cast<std::size_t>(i)])
                next.push_back(faces[static_cast<std::size_t>(i)]);
        for (const Face& f : horizon) next.push_back(f);
        faces.swap(next);
    }

    // simplicial polytope: every face is a triangle, so the count IS the answer.
    // and it checks out against Euler: faces.size() == 2*(distinct hull V) - 4.
    return static_cast<int>(faces.size());
}

}  // namespace p0051
