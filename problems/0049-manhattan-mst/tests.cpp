#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0049::manhattan_mst;
using p0049::Point;
using p0049::ReferenceSolver;

namespace {

std::int64_t mst(const std::vector<Point>& v) { return manhattan_mst(v); }

// a random point set: n in [lo_n, hi_n], each coord in [-range, range]. a small
// range is on purpose -- it forces ties, coincident points, and accidental
// collinearity, which is exactly where the octant sweep and the (x-y)
// compression are easiest to get wrong.
std::vector<Point> random_set(cpf::Rng& rng, int lo_n, int hi_n,
                              std::int64_t range) {
    int n = static_cast<int>(rng.in_range(lo_n, hi_n));
    std::vector<Point> v(static_cast<std::size_t>(n));
    for (auto& p : v) {
        p.x = rng.in_range(-range, range);
        p.y = rng.in_range(-range, range);
    }
    return v;
}

// everything on one axis-aligned or diagonal line -- the degenerate case where
// the whole MST is a single chain.
std::vector<Point> collinear_set(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 40));
    std::vector<Point> v(static_cast<std::size_t>(n));
    int mode = static_cast<int>(rng.in_range(0, 2));
    for (auto& p : v) {
        std::int64_t t = rng.in_range(-30, 30);
        if (mode == 0) {
            p = {t, 7};  // horizontal
        } else if (mode == 1) {
            p = {-4, t};  // vertical
        } else {
            p = {t, t};  // 45-degree diagonal -- octant boundaries everywhere
        }
    }
    return v;
}

}  // namespace

int main() {
    // -- hand examples --

    // two points: the answer is their distance, nothing else.
    CPF_EQ(mst({{0, 0}, {5, 3}}), 8);

    // unit square: three unit edges span it, the diagonals never get used.
    CPF_EQ(mst({{0, 0}, {1, 0}, {0, 1}, {1, 1}}), 3);

    // four points on a line: the MST is the chain, weight = the full span split
    // into consecutive gaps, 2 + 3 + 4.
    CPF_EQ(mst({{0, 0}, {2, 0}, {5, 0}, {9, 0}}), 9);

    // coincident points: distance zero, so the tree costs nothing.
    CPF_EQ(mst({{7, 7}, {7, 7}, {7, 7}}), 0);

    // -- edges --

    CPF_EQ(mst({}), 0);            // empty
    CPF_EQ(mst({{3, -9}}), 0);     // single point
    CPF_EQ(mst({{5, 5}, {5, 5}, {5, 5}, {5, 5}, {5, 5}}), 0);  // all identical

    // negative coordinates, mixed signs -- the sweep must not assume anything
    // about the octant living in the first quadrant.
    CPF_EQ(mst({{-3, -3}, {-3, 4}, {2, -3}}), 12);  // 7 + 5 chain
    CPF_EQ(mst({{-1000000000, -1000000000}, {1000000000, 1000000000}}),
           4000000000LL);  // full coordinate span, well past int32

    // -- randomized differential against the O(n^2) prim oracle --
    //
    // thousands of small sets across four regimes. any disagreement prints the
    // seed that reproduces it, and the seed is the whole repro.
    ReferenceSolver ref;
    auto oracle = [&](const std::vector<Point>& v) { return ref.solve(v); };

    int batches = 0;
    batches += cpf::differential(
        5000, 1, [](cpf::Rng& r) { return random_set(r, 1, 50, 15); }, mst,
        oracle);  // wide-ish: normal spread, occasional ties
    batches += cpf::differential(
        5000, 100000,
        [](cpf::Rng& r) { return random_set(r, 1, 40, 3); }, mst,
        oracle);  // tiny range: coincident points everywhere
    batches += cpf::differential(
        3000, 200000,
        [](cpf::Rng& r) { return random_set(r, 1, 30, 4); }, mst,
        oracle);  // grid-aligned 9x9 lattice
    batches += cpf::differential(3000, 300000,
                                 [](cpf::Rng& r) { return collinear_set(r); },
                                 mst, oracle);  // degenerate chains

    // differential returns true only when a whole batch matched.
    CPF_CHECK(batches == 4);
    std::printf("differential: 16000 small sets vs O(n^2) prim, %d/4 batches\n",
                batches);

    return cpf::report();
}
