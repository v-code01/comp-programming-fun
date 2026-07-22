#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0051::convex_hull_faces;
using p0051::convex_hull_faces_brute;
using p0051::P;
using p0051::detail::orient;

namespace {

int faces(const std::vector<P>& v) { return convex_hull_faces(v); }

// general position, stated as one predicate: every 4 points bound a real tetra,
// orient != 0. that single check catches BOTH forbidden degeneracies -- 4
// coplanar points give volume 0 directly, and 3 collinear points give volume 0
// against ANY fourth (a flat "plane" has no side). so this is exactly the
// input guarantee the hull leans on, and the generator below refuses anything
// that fails it.
bool is_general_position(const std::vector<P>& v) {
    const int n = static_cast<int>(v.size());
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            for (int c = b + 1; c < n; ++c)
                for (int d = c + 1; d < n; ++d)
                    if (orient(v[static_cast<std::size_t>(a)],
                               v[static_cast<std::size_t>(b)],
                               v[static_cast<std::size_t>(c)],
                               v[static_cast<std::size_t>(d)]) == 0)
                        return false;
    return true;
}

}  // namespace

int main() {
    // -- hand examples --

    // a tetrahedron. four points, four triangular faces. Euler: 2*4 - 4 = 4.
    CPF_EQ(faces({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}), 4);

    // a triangular bipyramid: a base triangle in z=0 with an apex above and one
    // below. neither apex, nor the base, is coplanar with the other four. each
    // apex caps the three base edges -> 3 + 3 = 6 faces. the base triangle is
    // interior, not a face. Euler: 2*5 - 4 = 6.
    CPF_EQ(faces({{0, 0, 0}, {6, 0, 0}, {0, 6, 0}, {1, 1, 5}, {2, 2, -5}}), 6);

    // interior point. a tetra plus one point strictly inside it -- the inside
    // point never reaches the surface, so the hull is still the tetra: 4 faces.
    CPF_EQ(faces({{0, 0, 0}, {6, 0, 0}, {0, 6, 0}, {0, 0, 6}, {1, 1, 1}}), 4);

    // -- edges --

    // exactly four points is the smallest hull -- one insertion loop never runs.
    CPF_EQ(faces({{-3, -3, -3}, {3, -1, 0}, {0, 4, -2}, {1, 0, 5}}), 4);

    // two interior points, not one. hull stays the tetra.
    CPF_EQ(faces({{0, 0, 0},
                  {9, 0, 0},
                  {0, 9, 0},
                  {0, 0, 9},
                  {1, 1, 1},
                  {2, 1, 2}}),
           4);

    // the generator's rejection actually fires: four points on z=0 are coplanar,
    // so orient == 0 and is_general_position says no. this is the exact gate that
    // keeps the differential fed only general-position sets.
    CPF_CHECK(!is_general_position(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}));
    // three collinear points are caught by the same test -- volume 0 against a
    // fourth off the line.
    CPF_CHECK(!is_general_position(
        {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {5, 0, 3}}));
    // a real tetra passes.
    CPF_CHECK(is_general_position(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));

    // -- randomized differential vs the O(n^4) brute oracle --
    //
    // thousands of small general-position sets. the generator draws random
    // integer points and REJECTS any set with 4 coplanar points (is_general_
    // position == false) before it ever reaches the hull -- so every set the
    // hull sees is a legal input. small coordinate ranges are on purpose: they
    // force coplanar rejections and pack interior points, the two spots the
    // horizon walk and the tetra seeding are easiest to get wrong.
    cpf::Rng rng(918273645ULL);
    int sets = 0, diffs = 0, rejects = 0;
    int first_bad = -1;

    struct Regime {
        int lo_n, hi_n;
        std::int64_t range;
        int want;  // how many accepted sets to run in this regime
    };
    const Regime regimes[] = {
        {4, 8, 4, 3000},    // tight range: heavy coplanar rejection, interior pts
        {4, 12, 12, 3000},  // mid: bigger hulls, occasional rejection
        {5, 15, 60, 3000},  // wide: rejections rare, hulls with many faces
    };

    for (const Regime& reg : regimes) {
        int accepted = 0;
        while (accepted < reg.want) {
            const int n = static_cast<int>(rng.in_range(reg.lo_n, reg.hi_n));
            std::vector<P> pts(static_cast<std::size_t>(n));
            for (auto& q : pts) {
                q.x = rng.in_range(-reg.range, reg.range);
                q.y = rng.in_range(-reg.range, reg.range);
                q.z = rng.in_range(-reg.range, reg.range);
            }
            if (!is_general_position(pts)) {
                ++rejects;  // coplanar-4 (or collinear-3) -- not a legal input
                continue;
            }
            ++accepted;
            ++sets;

            const int got = convex_hull_faces(pts);
            const int want = convex_hull_faces_brute(pts);
            ++cpf::state().total;
            if (got != want) {
                ++cpf::state().failed;
                if (diffs == 0) first_bad = static_cast<int>(pts.size());
                ++diffs;
            }
        }
    }

    CPF_CHECK(diffs == 0);
    if (diffs) std::printf("  first diff on a set of size %d\n", first_bad);

    // the generator has to actually reject coplanar sets, or the "general
    // position" claim is untested. with these small ranges it fires plenty.
    CPF_CHECK(rejects > 0);

    std::printf(
        "differential: %d general-position sets vs O(n^4) brute, %d diffs; "
        "generator rejected %d coplanar-4 sets\n",
        sets, diffs, rejects);

    return cpf::report();
}
