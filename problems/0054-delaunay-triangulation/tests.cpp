#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0054::P;

namespace {

using i128 = __int128;

// exact in-circle SIGN (the solution's predicate only exposes the > 0 bool; the
// general-position check needs to catch the == 0 concyclic case too). assumes
// (a,b,c) is CCW. +1 d strictly inside, -1 strictly outside, 0 ON the circle.
int in_circle_sign(const P& a, const P& b, const P& c, const P& d) {
    const i128 ax = a.x - d.x, ay = a.y - d.y;
    const i128 bx = b.x - d.x, by = b.y - d.y;
    const i128 cx = c.x - d.x, cy = c.y - d.y;
    const i128 az = ax * ax + ay * ay;
    const i128 bz = bx * bx + by * by;
    const i128 cz = cx * cx + cy * cy;
    const i128 det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) +
                     az * (bx * cy - by * cx);
    return (det > 0) - (det < 0);
}

int orient_sign(const P& a, const P& b, const P& c) {
    const i128 v =
        (i128)(b.x - a.x) * (c.y - a.y) - (i128)(b.y - a.y) * (c.x - a.x);
    return (v > 0) - (v < 0);
}

// general position: no 3 collinear, no 4 concyclic. the solution's points-at-
// infinity predicates assume every test is strictly nonzero, so a degenerate
// input isn't "wrong output", it's outside the contract -- the generator must
// never emit one, and the oracle differential would be meaningless on one.
bool is_general_position(const std::vector<P>& pts) {
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j + 1; k < n; ++k) {
                const int o = orient_sign(pts[static_cast<std::size_t>(i)],
                                          pts[static_cast<std::size_t>(j)],
                                          pts[static_cast<std::size_t>(k)]);
                if (o == 0) return false;  // 3 collinear
                int a = i, b = j, c = k;
                if (o < 0) std::swap(b, c);  // CCW for the in-circle sign
                for (int m = 0; m < n; ++m) {
                    if (m == i || m == j || m == k) continue;
                    if (in_circle_sign(pts[static_cast<std::size_t>(a)],
                                       pts[static_cast<std::size_t>(b)],
                                       pts[static_cast<std::size_t>(c)],
                                       pts[static_cast<std::size_t>(m)]) == 0)
                        return false;  // 4 concyclic
                }
            }
    return true;
}

// running count of draws the generator threw back -- proof it actually rejects
// degeneracy rather than hoping random floats miss it.
long g_rejects = 0;

std::vector<P> random_gp_pointset(cpf::Rng& rng) {
    const int n = static_cast<int>(rng.in_range(3, 9));
    for (;;) {
        std::vector<P> pts(static_cast<std::size_t>(n));
        for (auto& p : pts)
            p = {rng.in_range(-60, 60), rng.in_range(-60, 60)};
        if (is_general_position(pts)) return pts;
        ++g_rejects;
    }
}

std::vector<std::pair<int, int>> tri(int a, int b, int c) {
    std::vector<std::pair<int, int>> e = {{a, b}, {a, c}, {b, c}};
    return e;
}

}  // namespace

int main() {
    // -- hand case: a single triangle. all 3 sides are Delaunay edges. --
    {
        std::vector<P> pts = {{0, 0}, {4, 0}, {0, 3}};
        CPF_EQ(p0054::delaunay_edge_sq_sum(pts),
               static_cast<std::int64_t>(50));  // 16 + 25 + 9
        CPF_EQ(p0054::delaunay_edges(pts), tri(0, 1, 2));
    }

    // -- hand case: a non-concyclic quad. Delaunay picks the shorter diagonal
    // over the longer one; the solution must agree with the oracle exactly. --
    {
        std::vector<P> pts = {{0, 0}, {10, 0}, {11, 4}, {0, 3}};
        CPF_CHECK(is_general_position(pts));
        CPF_EQ(p0054::delaunay_edges(pts), p0054::delaunay_edges_brute(pts));
        CPF_EQ(p0054::delaunay_edge_sq_sum(pts),
               p0054::delaunay_edge_sq_sum_brute(pts));
    }

    // -- hand case: an interior point (fans into the surrounding triangle). --
    {
        std::vector<P> pts = {{0, 0}, {20, 0}, {10, 18}, {10, 6}};
        CPF_CHECK(is_general_position(pts));
        CPF_EQ(p0054::delaunay_edges(pts), p0054::delaunay_edges_brute(pts));
    }

    // -- the general-position gate itself: it must REJECT degeneracy. --
    {
        std::vector<P> collinear = {{0, 0}, {1, 1}, {2, 2}};
        CPF_CHECK(!is_general_position(collinear));         // 3 in a line
        std::vector<P> square = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
        CPF_CHECK(!is_general_position(square));            // 4 on one circle
        std::vector<P> good = {{0, 0}, {5, 0}, {1, 4}, {8, 7}};
        CPF_CHECK(is_general_position(good));               // real quad, accepted
    }

    // -- the differential: solution edge set == oracle edge set, over thousands
    // of random general-position clouds. compares the full SET, not just the
    // summed value, so a single flipped diagonal is caught. --
    std::printf("differential: solution vs O(n^4) empty-circle oracle\n");
    cpf::differential(
        15000, 1, [](cpf::Rng& rng) { return random_gp_pointset(rng); },
        [](const std::vector<P>& pts) { return p0054::delaunay_edges(pts); },
        [](const std::vector<P>& pts) {
            return p0054::delaunay_edges_brute(pts);
        });
    std::printf("  generator rejected %ld non-general-position draws\n",
                g_rejects);

    return cpf::report();
}
