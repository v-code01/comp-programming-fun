#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <numeric>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0042::Pt;
using ll = std::int64_t;

ll fast(const std::vector<Pt>& a, const std::vector<Pt>& b) {
    return p0042::solve(a, b);
}
ll brute(const std::vector<Pt>& a, const std::vector<Pt>& b) {
    return p0042::solve_ref(a, b);
}

// a COMPATIBLE configuration, by construction: x nondecreasing, y nonincreasing,
// in both lists. build it by sorting the coordinates independently and zipping --
// then the ordering hypothesis can't be violated by accident, and the differential
// is testing the determinant, not the generator.
//
// ties are left in on purpose. equal starts (or equal ends) mean two identical
// rows (or columns), so the determinant is 0 -- and the brute walk agrees, since
// two paths sharing an endpoint share a lattice point. that agreement is worth
// testing, not worth generating around.
std::vector<Pt> sorted_points(cpf::Rng& rng, int k, int lo, int hi) {
    std::vector<ll> xs = rng.vec(static_cast<std::size_t>(k), lo, hi);
    std::vector<ll> ys = rng.vec(static_cast<std::size_t>(k), lo, hi);
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end(), std::greater<ll>());
    std::vector<Pt> p(static_cast<std::size_t>(k));
    for (int i = 0; i < k; ++i)
        p[static_cast<std::size_t>(i)] = Pt{static_cast<int>(xs[static_cast<std::size_t>(i)]),
                                            static_cast<int>(ys[static_cast<std::size_t>(i)])};
    return p;
}

// the classic pin. starts (0,1),(1,0), ends (n-1,n),(n,n-1) -- the two-path ballot
// setup, the one every LGV writeup opens with. write m = n-1. the matrix is
//     [[C(2m,m), C(2m,m+1)], [C(2m,m-1), C(2m,m)]]
// so the determinant is C(2m,m)^2 - C(2m,m+1)*C(2m,m-1), and the vandermonde-ish
// collapse of that is
//     (2m+1) * Catalan(m)^2 .
// for n = 1..6: 1, 3, 20, 175, 1764, 19404. an independent number to hit, with no
// determinant anywhere near it.
ll catalan(int m) {
    ll c = 1;  // C(2m,m)/(m+1), built with exact integer arithmetic -- m stays small.
    for (int i = 0; i < m; ++i) c = c * (2 * m - i) / (i + 1);
    return c / (m + 1);
}

}  // namespace

int main() {
    using namespace p0042;

    // ---- k = 1. the determinant of a 1x1 matrix is the binomial, nothing else. ----

    CPF_EQ(fast({{0, 0}}, {{3, 2}}), 10);   // C(5,3) = 10.
    CPF_EQ(brute({{0, 0}}, {{3, 2}}), 10);  // and the walk agrees.
    CPF_EQ(fast({{2, 7}}, {{2, 7}}), 1);    // start == end. one path: stand still.
    CPF_EQ(fast({{5, 5}}, {{1, 1}}), 0);    // end behind the start. no monotone path.
    CPF_EQ(fast({{5, 5}}, {{9, 1}}), 0);    // end below the start. still nothing.

    // the far corner of the table. C(n,0) = 1: a straight line right, or straight
    // up, is the only monotone path there is. this is the check that the factorial
    // tables actually reach 4e5 -- undersize them and it reads off the end.
    CPF_EQ(fast({{0, 0}}, {{200000, 0}}), 1);
    CPF_EQ(fast({{0, 0}}, {{0, 200000}}), 1);
    // C(n,r) == C(n,n-r), so a 2e5 x 2e5 corner and its transpose agree -- and the
    // count is nonzero, which an off-by-one in the inverse table would break.
    CPF_EQ(fast({{0, 0}}, {{200000, 199999}}), fast({{0, 0}}, {{199999, 200000}}));
    CPF_CHECK(fast({{0, 0}}, {{200000, 200000}}) != 0);

    // ---- k = 2, two paths that have room. ----

    // starts (0,1),(1,0), ends (1,2),(2,1). M = [[2,1],[1,2]], det = 3.
    {
        const std::vector<Pt> a = {{0, 1}, {1, 0}};
        const std::vector<Pt> b = {{1, 2}, {2, 1}};
        CPF_EQ(fast(a, b), 3);
        CPF_EQ(brute(a, b), 3);
    }

    // a corridor that costs you paths but doesn't close: a_2 sits on two of the
    // three routes a_1 could take, so 3 families become 2.
    {
        const std::vector<Pt> a = {{0, 0}, {1, 0}};
        const std::vector<Pt> b = {{1, 2}, {2, 0}};
        CPF_EQ(fast(a, b), 2);
        CPF_EQ(brute(a, b), 2);
    }

    // ---- forced crossings -> 0. ----

    // everything on one row. path 1 is (0,0)->(1,0), path 2 is (1,0)->(2,0). the
    // only route for each is the straight one, and they collide at (1,0). every
    // entry of M is 1, and det [[1,1],[1,1]] = 0. the zero is the crossing.
    {
        const std::vector<Pt> a = {{0, 0}, {1, 0}};
        const std::vector<Pt> b = {{1, 0}, {2, 0}};
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // three on a row. same story, 3x3.
    {
        const std::vector<Pt> a = {{0, 0}, {1, 0}, {2, 0}};
        const std::vector<Pt> b = {{2, 0}, {3, 0}, {4, 0}};
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // duplicate starts -- two identical rows. the paths share a_1 = a_2 before
    // they've taken a single step.
    {
        const std::vector<Pt> a = {{1, 1}, {1, 1}};
        const std::vector<Pt> b = {{2, 3}, {3, 2}};
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // duplicate ends -- two identical columns.
    {
        const std::vector<Pt> a = {{0, 2}, {2, 0}};
        const std::vector<Pt> b = {{3, 3}, {3, 3}};
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // one end is unreachable from its start -- b_2 is strictly below a_2. that row
    // of M is not all-zero (b_1 is still reachable), so the zero has to come out of
    // the elimination, not out of an early bail.
    {
        const std::vector<Pt> a = {{0, 3}, {2, 1}};
        const std::vector<Pt> b = {{1, 4}, {3, 0}};  // a_2=(2,1) -> b_2=(3,0): y drops.
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // nothing reachable at all -- both ends behind both starts. M is the zero
    // matrix, no pivot in column 0.
    {
        const std::vector<Pt> a = {{5, 6}, {6, 5}};
        const std::vector<Pt> b = {{0, 1}, {1, 0}};
        CPF_EQ(fast(a, b), 0);
        CPF_EQ(brute(a, b), 0);
    }

    // ---- the ballot pin: with m = n-1, the count is (2m+1) * Catalan(m)^2. ----
    {
        for (int n = 1; n <= 6; ++n) {
            const int m = n - 1;
            const std::vector<Pt> a = {{0, 1}, {1, 0}};
            const std::vector<Pt> b = {{n - 1, n}, {n, n - 1}};
            const ll want = (2 * m + 1) * catalan(m) * catalan(m);
            CPF_EQ(fast(a, b), want);
            CPF_EQ(brute(a, b), want);  // coords stay inside the oracle's 0..7 grid.
        }
    }

    // ---- the compatibility condition, tested and not assumed. ----
    //
    // LGV gives a SIGNED sum. it only equals the plain count because no non-identity
    // pairing admits a disjoint family. so ask the brute walk directly: over random
    // compatible configs, for every sigma != id, the disjoint-family count must be 0.
    // if this ever fires, the ordering hypothesis is wrong and the determinant is
    // counting something else.
    {
        int checked = 0, violations = 0;
        for (int k = 2; k <= 3; ++k) {
            cpf::Rng rng(0x0042C0 + static_cast<std::uint64_t>(k));
            for (int c = 0; c < 25000; ++c) {
                const auto a = sorted_points(rng, k, 0, 5);
                const auto b = sorted_points(rng, k, 0, 5);
                std::vector<int> perm(static_cast<std::size_t>(k));
                std::iota(perm.begin(), perm.end(), 0);
                while (std::next_permutation(perm.begin(), perm.end())) {
                    ++checked;
                    if (solve_ref(a, b, perm) != 0) ++violations;
                }
            }
        }
        CPF_CHECK(violations == 0);
        std::printf("compatibility: %d non-identity pairings, %d admit a disjoint family\n",
                    checked, violations);
    }

    // ---- randomized differential vs the brute path enumeration. ----

    int cases = 0, zeros = 0, nonzeros = 0, diffs = 0, first_bad = -1;

    auto sweep = [&](std::uint64_t seed, int n, int k, int lo, int hi) {
        cpf::Rng rng(seed);
        for (int c = 0; c < n; ++c) {
            const auto a = sorted_points(rng, k, lo, hi);
            const auto b = sorted_points(rng, k, lo, hi);
            const ll want = brute(a, b);
            CPF_CHECK(want >= 0);  // -1 means the oracle was handed coords it can't pack.
            const ll got = fast(a, b);
            ++cases;
            if (want == 0) ++zeros;
            else ++nonzeros;
            if (got != want) {
                if (diffs == 0) first_bad = c;
                ++diffs;
            }
        }
    };

    // k=1 -- the pure binomial, plus every unreachable orientation.
    sweep(0x00421, 20000, 1, 0, 7);
    // k=2 on a tight grid -- crossings everywhere, so plenty of honest zeros.
    sweep(0x00422, 60000, 2, 0, 4);
    // k=2 with more room -- longer paths, bigger counts.
    sweep(0x00423, 30000, 2, 0, 7);
    // k=3 -- three-way collisions, where a sign slip in the 3x3 elimination shows.
    sweep(0x00424, 50000, 3, 0, 4);
    // k=3 with room to breathe -- the nonzero counts get big enough to be real.
    sweep(0x00425, 20000, 3, 0, 6);
    // k=3 packed into a 4x4 -- almost everything here is forced to cross.
    sweep(0x00426, 20000, 3, 0, 3);

    CPF_CHECK(diffs == 0);
    if (diffs) std::printf("  first diff in a sweep at case %d\n", first_bad);
    std::printf("differential vs brute enumeration: %d cases, %d zeros, %d nonzero, %d diffs\n",
                cases, zeros, nonzeros, diffs);

    return cpf::report();
}
