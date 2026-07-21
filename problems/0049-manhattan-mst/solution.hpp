#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace p0049 {

// manhattan minimum spanning tree.
//
// n points in the plane. the edge between two of them weighs their manhattan
// distance |dx| + |dy|. we want the total weight of the MST over the complete
// graph -- all C(n,2) edges.
//
// the naive thing builds every edge: O(n^2) to enumerate, O(n^2 log n) to
// kruskal. at n = 1e5 that's 5e9 edges. dead on arrival.
//
// lower bound. reading the points is Omega(n). the sweep sorts, and sorting in
// the comparison model is Omega(n log n). so O(n log n) is the floor and we hit
// it -- the candidate reduction below drops the edge count from n^2 to <= 4n.
//
// -- the candidate theorem --
// slice the directions around a point into eight 45-degree octants. the claim:
// for any point p and any octant, at most one edge out of p in that octant can
// belong to the MST -- and it's the edge to p's NEAREST point in that octant.
//
// why. take two points s and t in the same 45-degree octant of a third point p,
// with t the nearer of the two (dist(p,t) <= dist(p,s)). inside one octant the
// L1 geometry forces dist(t,s) <= dist(p,s). so on the triangle {p,t,s} the edge
// p-s is a maximum. the MST cycle property throws out a max-weight edge on any
// cycle -- so p-s is never forced into the MST. only the nearest survives.
// (guibas & stolfi. the 45-degree width is load-bearing: a wider cone breaks the
// dist(t,s) <= dist(p,s) step.)
//
// so the true MST is a subset of { (p, nearest-in-octant(p)) }. that's <= 8n
// edges. kruskal over them equals the full MST -- same weight, because an MST
// lives inside this set.
//
// -- four sweeps, not eight --
// every edge is undirected. orient it by (x+y): the endpoint with the smaller
// x+y "looks toward" the larger. that direction always lands in the closed half
// [-45, 135] degrees -- half the octants. so we only sweep the four 45-degree
// octants that tile that half:  [-45,0]  [0,45]  [45,90]  [90,135].  every MST
// edge, oriented low-to-high, sits in exactly one of them.
//
// -- one sweep --
// take the base octant [0,45]:  dx >= dy >= 0. it splits into two independent
// half-plane conditions plus a linear cost:
//   dy >= 0     <=>   y_q >= y_p            -- the "above" half
//   dx >= dy    <=>   (x-y)_q >= (x-y)_p    -- the "shallow enough" half
//   dist(p,q) = dx + dy = (x+y)_q - (x+y)_p -- so nearest = min (x+y)_q.
// (careful: y_q>=y_p and (x-y)_q>=(x-y)_p are BOTH needed. you cannot fold them
// into one x+y test -- subtracting two inequalities is not a valid step, and a
// point with dy<0 would sneak into the octant and win the query. that exact
// mistake makes the MST come out too heavy.)
//
// so: process points by DESCENDING y. when we reach p, every inserted q already
// has y_q >= y_p -- the "above" half comes for free. among those we still need
// (x-y)_q >= (x-y)_p and the smallest (x+y)_q -- one suffix-min query on a
// fenwick keyed by compressed (x-y), holding min (x+y). insert p, move on.
// O(n log n) per sweep. ties in y break by descending x so a purely horizontal
// neighbour (dy=0, dx>0) is already inserted when its partner asks for it.
//
// -- the four transforms --
// each sweep runs the SAME base-octant code on transformed coordinates. we pick
// T so the base octant [0,45] in T-space pulls back to the octant we want:
//   s=0  (x,  y)   -> [0,45]
//   s=1  (y,  x)   -> [45,90]    reflect across y=x  (theta -> 90-theta)
//   s=2  (y, -x)   -> [90,135]   rotate -90          (theta -> theta-90)
//   s=3  (x, -y)   -> [-45,0]    reflect across x    (theta -> -theta)
// all four are L1-isometries, so weights survive untouched. their images tile
// [-45,135]. done.
//
// -- overflow --
// coords in [-1e9, 1e9]. x+y and x-y reach 2e9 -- past int32, so everything is
// int64. an edge weight reaches 4e9; the tree sums <= (n-1) of them, ~4e14,
// still int64. no int anywhere on the value path.

struct Point {
    std::int64_t x;
    std::int64_t y;
};

namespace detail {

inline std::int64_t iabs(std::int64_t v) { return v < 0 ? -v : v; }

inline std::int64_t manhattan(const Point& a, const Point& b) {
    return iabs(a.x - b.x) + iabs(a.y - b.y);
}

// fenwick over compressed keys, holding the prefix-min of (value, point-index).
// 1-indexed. update pins a value at one key; query returns the min over a prefix
// [1..r] and the index that owns it. we compress (x-y) so that a LARGER (x-y)
// gets a SMALLER rank -- then the octant's "(x-y)_q >= (x-y)_p" suffix turns
// into the prefix [1..rank(p)] the fenwick answers natively.
struct FenwickMin {
    static constexpr std::int64_t kInf = std::numeric_limits<std::int64_t>::max();
    std::vector<std::int64_t> val;
    std::vector<int> idx;

    explicit FenwickMin(int n) : val(n + 1, kInf), idx(n + 1, -1) {}

    void update(int i, std::int64_t v, int id) {
        for (; i < static_cast<int>(val.size()); i += i & -i) {
            if (v < val[i]) {
                val[i] = v;
                idx[i] = id;
            }
        }
    }

    // min value over [1..r] and its owner index (-1 if the prefix is empty).
    std::pair<std::int64_t, int> query(int r) const {
        std::int64_t best = kInf;
        int who = -1;
        for (; r > 0; r -= r & -r) {
            if (val[r] < best) {
                best = val[r];
                who = idx[r];
            }
        }
        return {best, who};
    }
};

struct Edge {
    std::int64_t w;
    int u;
    int v;
};

// union-find with path compression + union by size. near-constant per op, so
// kruskal over ~4n edges is dominated by the sort, not the merges.
struct Dsu {
    std::vector<int> parent;
    std::vector<int> size;

    explicit Dsu(int n) : parent(n), size(n, 1) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (size[a] < size[b]) std::swap(a, b);
        parent[b] = a;
        size[a] += size[b];
        return true;
    }
};

// one octant sweep. tx/ty are the transformed coordinates; orig feeds the true
// edge weight. appends one candidate edge per point that has a neighbour in the
// base octant [0,45] of the transformed frame.
inline void sweep_octant(const std::vector<Point>& orig,
                         const std::vector<std::int64_t>& tx,
                         const std::vector<std::int64_t>& ty,
                         std::vector<Edge>& edges) {
    const int n = static_cast<int>(orig.size());

    // compress the keys (x-y), largest first so bigger key -> smaller rank.
    std::vector<std::int64_t> key(n);
    for (int i = 0; i < n; ++i) key[i] = tx[i] - ty[i];
    std::vector<std::int64_t> uniq = key;
    std::sort(uniq.begin(), uniq.end(), std::greater<std::int64_t>());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    auto rank_of = [&](std::int64_t k) {
        // first position with element <= k in the descending array -- since k is
        // present, that's k's own slot. +1 for 1-indexing.
        int pos = static_cast<int>(
            std::lower_bound(uniq.begin(), uniq.end(), k,
                             std::greater<std::int64_t>()) -
            uniq.begin());
        return pos + 1;
    };

    // process by descending y: everything already inserted has y_q >= y_p, which
    // is the octant's "dy >= 0" half. ties break by descending x so a horizontal
    // neighbour is in the fenwick before its partner queries.
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return ty[a] != ty[b] ? ty[a] > ty[b] : tx[a] > tx[b];
    });

    FenwickMin fen(static_cast<int>(uniq.size()));
    for (int p : order) {
        int r = rank_of(key[p]);
        auto [best, who] = fen.query(r);  // nearest ahead-and-in-octant, or -1
        if (who != -1) {
            edges.push_back({manhattan(orig[p], orig[who]), p, who});
        }
        fen.update(r, tx[p] + ty[p], p);
    }
}

}  // namespace detail

// total weight of the manhattan MST. 0 for n <= 1. O(n log n).
inline std::int64_t manhattan_mst(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    if (n <= 1) return 0;

    std::vector<detail::Edge> edges;
    edges.reserve(static_cast<std::size_t>(4) * n);

    std::vector<std::int64_t> tx(n), ty(n);
    for (int s = 0; s < 4; ++s) {
        for (int i = 0; i < n; ++i) {
            const std::int64_t x = pts[i].x, y = pts[i].y;
            switch (s) {
                case 0: tx[i] = x;  ty[i] = y;  break;  // [0,45]
                case 1: tx[i] = y;  ty[i] = x;  break;  // [45,90]
                case 2: tx[i] = y;  ty[i] = -x; break;  // [90,135]
                default: tx[i] = x; ty[i] = -y; break;  // [-45,0]
            }
        }
        detail::sweep_octant(pts, tx, ty, edges);
    }

    // kruskal over the <= 4n candidates. the sort is the whole cost.
    std::sort(edges.begin(), edges.end(),
              [](const detail::Edge& a, const detail::Edge& b) {
                  return a.w < b.w;
              });

    detail::Dsu dsu(n);
    std::int64_t total = 0;
    int used = 0;
    for (const auto& e : edges) {
        if (dsu.unite(e.u, e.v)) {
            total += e.w;
            if (++used == n - 1) break;  // spanning tree complete.
        }
    }
    return total;
}

}  // namespace p0049
