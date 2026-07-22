#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0054 {

// Delaunay triangulation -- sum of squared edge lengths.
//
// n points, integer coords, general position: no 3 collinear, no 4 concyclic.
// build the Delaunay triangulation, take its edge set, sum (dx^2 + dy^2) over
// every edge. one exact integer, fits int64.
//
// -- the size of the answer, and the floor --
// reading the points is Omega(n). a triangulation of n points with h of them on
// the convex hull has, by Euler (V - E + F = 2, counting the outer face and
// 2E = 3(F-1) + h since every bounded face is a triangle and the outer face has
// h edges):
//   triangles = 2n - 2 - h,   edges = 3n - 3 - h.
// so the output is Theta(n) edges -- linear in the input, and the reading floor
// is tight. the brute oracle that defines "Delaunay edge" from scratch is O(n^4)
// (every pair, every witness third point, every other point). we don't ship that
// -- we ship Bowyer-Watson and diff against it.
//
// -- Bowyer-Watson, incremental, O(n^2) --
// wrap the whole point set in a super-triangle, then insert points one at a time.
// for a new point p, a triangle is BAD if p is strictly inside its circumcircle.
// the bad triangles form one star-shaped cavity around p. delete them, keep the
// cavity boundary -- each boundary edge belongs to exactly one bad triangle --
// and fan a fresh triangle from p to every boundary edge. at the end drop every
// triangle that still touches a super vertex; the survivors are the Delaunay
// triangulation of the real points.
// each insertion scans every current triangle for badness: O(F) = O(n) triangles
// (F = 2n-2-h). n insertions -> O(n^2), ~8e6 in-circle calls at n=2000. a spatial
// walk to the containing triangle would cut the scan to O(n log n) expected; we
// keep the flat scan on purpose -- it is easier to prove correct, and exactness
// is the thing this problem actually rewards.
//
// -- the in-circle predicate, exact, and the magnitude bound --
// for a CCW triangle (a,b,c), point d is strictly inside its circumcircle iff
//   det | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
//       | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |  > 0.
//       | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
// coords are <= 1e4, so a difference is <= 2e4, a squared term <= 8e8, and each
// of the six determinant products is <= (2e4)(2e4)(8e8) ~ 3.2e17 -- the sum is
// <= ~2e21. that overflows int64 (9.2e18) but sits far inside __int128 (1.7e38).
// no float anywhere: the sign has to be exact or the triangulation flips an edge.
//
// -- the super-triangle, handled as points at INFINITY --
// here is the trap. a near-collinear integer triangle in [-1e4,1e4]^2 can have a
// circumradius up to ~1e13 (min area 1/2, sides up to ~2.8e4). so NO finite
// super-triangle both encloses every such circumcircle AND keeps the determinant
// inside int128 -- to dominate a 1e13 radius you'd push super coords past 1e13,
// and the determinant would then need ~1e52, past int128. a modest finite
// super-triangle would silently give the wrong answer on thin near-hull
// triangles, exactly the ones a 2000-point cloud is full of.
//
// so the three super vertices are treated as points at INFINITY in three spread
// directions. their finite coords below are only DIRECTIONS -- the magnitude is
// symbolic (T -> infinity). taking that limit, the in-circle determinant of a
// super-incident triangle collapses to a lower-degree test:
//   * one super vertex (triangle a,b,SUPER with a,b real): the circle through two
//     fixed points and a third going to infinity becomes the LINE a-b, and its
//     interior becomes the half-plane on the super vertex's side. so p is inside
//     iff p is strictly left of the real edge a->b: orient(a,b,p) > 0.
//   * two super vertices (triangle a,SUPER1,SUPER2, a real): the leading term of
//     the determinant is  T^3 * cross(a - p, V),  where
//     V = e1*|e2|^2 - e2*|e1|^2  from the two super directions e1,e2 in the
//     triangle's cyclic order. so p is inside iff cross(a - p, V) > 0.
//   * three super vertices (only the initial triangle, before any real insertion):
//     every real point is inside by construction -- return true, so the first
//     insertion splits it into three.
// these reductions are the exact T->infinity limits (checked numerically against
// the finite determinant: the disagreement count falls ~10x per 10x T, pure
// finite-T boundary effect). DT of {real points + 3 points at infinity}, with the
// infinite-incident triangles stripped, is exactly DT of the real points -- hull
// edges and all. no super-triangle-size tuning, no lost thin triangles.
//
// -- why the super directions are BIG integers --
// the two-super test cross(a - p, V) could tie (== 0) if a - p were parallel to
// V. we pick super directions of magnitude ~6e4, which makes |V| ~ 1e14. the
// shortest integer vector parallel to V is then far longer than any a - p (whose
// length is <= ~2.8e4), so a - p is parallel to V only when a == p -- which never
// happens for two distinct points. the one-super test can't tie either: it needs
// three real points collinear, which general position forbids. so every
// predicate is strictly nonzero on legal input -- no perturbation, no tie-break.

struct P {
    std::int64_t x, y;
};

namespace detail {

using i128 = __int128;

// the three super vertices, as DIRECTIONS to infinity (CCW, well spread: up,
// lower-left, lower-right). magnitude ~6e4 on purpose -- see the note above on
// why big directions kill the two-super tie. these numbers never enter a squared
// determinant term; they appear only as directions and in the bounded V cross
// product, so int128 has room to spare.
inline constexpr P kSuper[3] = {{0, 60000}, {-60000, -20000}, {60000, -20000}};

// CCW check: orient((0,60000),(-60000,-20000),(60000,-20000)) = 9.6e9 > 0, so the
// fixed order 0,1,2 is already counter-clockwise -- no runtime reordering needed.

struct Tri {
    int a, b, c;  // CCW; an index >= n names a super vertex (index - n in kSuper)
};

// sign of orient(a,b,c) = sign of cross(b-a, c-a). +1 CCW, -1 CW, 0 collinear.
inline int orient(const P& a, const P& b, const P& c) {
    const i128 v = (i128)(b.x - a.x) * (c.y - a.y) -
                   (i128)(b.y - a.y) * (c.x - a.x);
    return (v > 0) - (v < 0);
}

// the all-real in-circle determinant. CCW triangle (a,b,c); d strictly inside
// circumcircle iff > 0. __int128 is load-bearing: the sum reaches ~2e21.
inline bool in_circle_real(const P& a, const P& b, const P& c, const P& d) {
    const i128 ax = a.x - d.x, ay = a.y - d.y;
    const i128 bx = b.x - d.x, by = b.y - d.y;
    const i128 cx = c.x - d.x, cy = c.y - d.y;
    const i128 az = ax * ax + ay * ay;
    const i128 bz = bx * bx + by * by;
    const i128 cz = cx * cx + cy * cy;
    const i128 det = ax * (by * cz - bz * cy)    //
                    - ay * (bx * cz - bz * cx)   //
                    + az * (bx * cy - by * cx);
    return det > 0;
}

}  // namespace detail

// build DT and return its edges as sorted, unique, normalized (u < v) index
// pairs. wants n >= 3 real points in general position -- the problem guarantees
// it. O(n^2).
inline std::vector<std::pair<int, int>> delaunay_edges(
    const std::vector<P>& pts) {
    using detail::in_circle_real;
    using detail::kSuper;
    using detail::orient;
    using detail::Tri;

    const int n = static_cast<int>(pts.size());
    std::vector<std::pair<int, int>> edges;
    if (n < 3) return edges;  // no triangle to speak of; keeps callers honest.

    // coord(i): real point for i < n, else the super DIRECTION for i - n.
    auto coord = [&](int i) -> P {
        return i < n ? pts[static_cast<std::size_t>(i)] : kSuper[i - n];
    };
    auto is_super = [&](int i) { return i >= n; };

    // p strictly inside the circumcircle of CCW triangle t. super vertices are at
    // infinity, so the test drops in degree with each super vertex present.
    auto bad = [&](const Tri& t, const P& p) -> bool {
        const int v[3] = {t.a, t.b, t.c};
        int sk = -1, sc = 0;
        for (int k = 0; k < 3; ++k)
            if (is_super(v[k])) {
                ++sc;
                sk = k;  // remembers a super slot (used when exactly one)
            }
        if (sc == 0) return in_circle_real(coord(v[0]), coord(v[1]),
                                           coord(v[2]), p);
        if (sc == 3) return true;  // initial triangle -- every real point inside
        if (sc == 1) {
            // rotate so the super vertex is last: reals are the two after it in
            // cyclic order. inside iff p is strictly left of that real edge.
            const P a = coord(v[(sk + 1) % 3]);
            const P b = coord(v[(sk + 2) % 3]);
            return orient(a, b, p) > 0;
        }
        // sc == 2: rotate so the lone real vertex is first, supers follow in
        // cyclic order. leading term is T^3 * cross(a - p, V).
        int rk = 0;
        for (int k = 0; k < 3; ++k)
            if (!is_super(v[k])) rk = k;
        const P a = coord(v[rk]);
        const P e1 = coord(v[(rk + 1) % 3]);
        const P e2 = coord(v[(rk + 2) % 3]);
        const detail::i128 e1n = (detail::i128)e1.x * e1.x + (detail::i128)e1.y * e1.y;
        const detail::i128 e2n = (detail::i128)e2.x * e2.x + (detail::i128)e2.y * e2.y;
        const detail::i128 vx = (detail::i128)e1.x * e2n - (detail::i128)e2.x * e1n;
        const detail::i128 vy = (detail::i128)e1.y * e2n - (detail::i128)e2.y * e1n;
        const detail::i128 adx = a.x - p.x, ady = a.y - p.y;
        const detail::i128 crs = adx * vy - ady * vx;  // cross(a - p, V)
        return crs > 0;
    };

    // seed the one super-triangle. CCW by the fixed order (checked above).
    std::vector<Tri> tris;
    tris.reserve(static_cast<std::size_t>(2) * n + 8);
    tris.push_back({n + 0, n + 1, n + 2});

    // horizon bookkeeping, same trick as the 3D hull: stamp[u*W + v] == epoch
    // means u->v is a directed edge of a bad triangle this round. W = n + 3
    // covers real and super indices. bump the epoch instead of clearing -- old
    // stamps go stale for free.
    const int W = n + 3;
    std::vector<int> stamp(static_cast<std::size_t>(W) * W, 0);
    int epoch = 0;

    std::vector<char> is_bad;
    std::vector<Tri> cavity;  // fresh triangles fanned from p this round
    std::vector<Tri> next;

    for (int p = 0; p < n; ++p) {
        ++epoch;
        const P pp = pts[static_cast<std::size_t>(p)];
        const int tc = static_cast<int>(tris.size());
        is_bad.assign(static_cast<std::size_t>(tc), 0);

        // mark bad triangles and stamp their directed edges.
        for (int i = 0; i < tc; ++i) {
            const Tri& t = tris[static_cast<std::size_t>(i)];
            if (bad(t, pp)) {
                is_bad[static_cast<std::size_t>(i)] = 1;
                stamp[static_cast<std::size_t>(t.a) * W + t.b] = epoch;
                stamp[static_cast<std::size_t>(t.b) * W + t.c] = epoch;
                stamp[static_cast<std::size_t>(t.c) * W + t.a] = epoch;
            }
        }

        // cavity boundary = a bad triangle's directed edge u->v whose reverse
        // v->u is NOT a bad edge (the triangle across it survives). fan p to each
        // as (u, v, p): p is left of the CCW boundary edge, so (u,v,p) stays CCW.
        cavity.clear();
        auto emit = [&](int u, int v) {
            if (stamp[static_cast<std::size_t>(v) * W + u] != epoch)
                cavity.push_back({u, v, p});
        };
        for (int i = 0; i < tc; ++i) {
            if (!is_bad[static_cast<std::size_t>(i)]) continue;
            const Tri& t = tris[static_cast<std::size_t>(i)];
            emit(t.a, t.b);
            emit(t.b, t.c);
            emit(t.c, t.a);
        }

        // drop the bad, keep the rest, append the fan.
        next.clear();
        next.reserve(static_cast<std::size_t>(tc) + cavity.size());
        for (int i = 0; i < tc; ++i)
            if (!is_bad[static_cast<std::size_t>(i)])
                next.push_back(tris[static_cast<std::size_t>(i)]);
        for (const Tri& t : cavity) next.push_back(t);
        tris.swap(next);
    }

    // strip every triangle still touching a super vertex; collect the real edges.
    for (const Tri& t : tris) {
        if (is_super(t.a) || is_super(t.b) || is_super(t.c)) continue;
        int e[3][2] = {{t.a, t.b}, {t.b, t.c}, {t.c, t.a}};
        for (auto& pr : e) {
            int u = pr[0], v = pr[1];
            if (u > v) std::swap(u, v);
            edges.emplace_back(u, v);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

// sum of squared lengths over the Delaunay edge set. exact int64.
inline std::int64_t delaunay_edge_sq_sum(const std::vector<P>& pts) {
    const auto edges = delaunay_edges(pts);
    std::int64_t sum = 0;
    for (const auto& e : edges) {
        const P& u = pts[static_cast<std::size_t>(e.first)];
        const P& v = pts[static_cast<std::size_t>(e.second)];
        const std::int64_t dx = u.x - v.x, dy = u.y - v.y;
        sum += dx * dx + dy * dy;  // <= 2*(2e4)^2 = 8e8 each, O(n) terms: no overflow
    }
    return sum;
}

}  // namespace p0054
