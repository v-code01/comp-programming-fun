#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0052::i64;
using p0052::kMod;

namespace {

std::vector<i64> rand_poly(int n, cpf::Rng& rng) {
    std::vector<i64> a(static_cast<std::size_t>(n));
    for (auto& x : a) x = rng.in_range(0, kMod - 1);
    return a;
}

std::vector<i64> rand_points(int m, cpf::Rng& rng) {
    std::vector<i64> p(static_cast<std::size_t>(m));
    for (auto& x : p) x = rng.in_range(0, kMod - 1);
    return p;
}

}  // namespace

int main() {
    using p0052::horner_multipoint;
    using p0052::multipoint_eval;

    // ---- hand examples worked out by hand ----
    {
        // A = 1 + 2x + 3x^2 at x = 2. Horner: 3, 3*2+2=8, 8*2+1=17. -> 17.
        CPF_EQ(multipoint_eval({1, 2, 3}, {2}), (std::vector<i64>{17}));

        // constant poly A = 7 -- evaluates to 7 at every point, x=0 included.
        {
            std::vector<i64> A{7};
            std::vector<i64> pts{0, 1, 5, 998244352, 123456};
            std::vector<i64> want(pts.size(), 7);
            CPF_EQ(multipoint_eval(A, pts), want);
        }

        // A = x = [0, 1] -- returns the points themselves.
        {
            std::vector<i64> A{0, 1};
            std::vector<i64> pts{0, 1, 2, 100, 998244352};
            CPF_EQ(multipoint_eval(A, pts), pts);
        }

        // A = [1,2,3] at several points, cross-checked against Horner.
        {
            std::vector<i64> A{1, 2, 3};
            std::vector<i64> pts{0, 1, 2, 3, 10};
            // by hand: 1, 6, 17, 34, 321.
            CPF_EQ(multipoint_eval(A, pts),
                   (std::vector<i64>{1, 6, 17, 34, 321}));
            CPF_EQ(multipoint_eval(A, pts), horner_multipoint(A, pts));
        }
    }

    // ---- edges ----
    {
        // n = 1, m = 1: constant a_0 at one point.
        CPF_EQ(multipoint_eval({42}, {5}), (std::vector<i64>{42}));

        // x = 0 pulls out a_0 exactly.
        {
            std::vector<i64> A{9, 8, 7, 6};
            CPF_EQ(multipoint_eval(A, {0}), (std::vector<i64>{9}));
        }

        // repeated points -- every copy must land the same value.
        {
            std::vector<i64> A{5, 4, 3, 2, 1};
            std::vector<i64> pts{7, 7, 7, 2, 2, 0, 0, 0};
            CPF_EQ(multipoint_eval(A, pts), horner_multipoint(A, pts));
        }

        // m < n: few points, tall polynomial (forces the real A-mod-root divide).
        {
            std::vector<i64> A{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};  // n = 11
            std::vector<i64> pts{1, 2, 3};                        // m = 3
            CPF_EQ(multipoint_eval(A, pts), horner_multipoint(A, pts));
        }

        // m > n: many points, short polynomial (root reduction is identity).
        {
            std::vector<i64> A{2, 3};  // n = 2
            std::vector<i64> pts{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
            CPF_EQ(multipoint_eval(A, pts), horner_multipoint(A, pts));
        }

        // coefficient near p and point near p -- overflow-sensitive corner.
        {
            std::vector<i64> A{kMod - 1, kMod - 2, 1, kMod - 5};
            std::vector<i64> pts{kMod - 1, kMod - 3, 0, 1};
            CPF_EQ(multipoint_eval(A, pts), horner_multipoint(A, pts));
        }
    }

    // ---- randomized differential vs Horner, thousands of small instances.
    //      n, m in 1..64, random coeffs and points, covering m<n, m>n, repeats,
    //      and x=0. any single mismatch is the tree or the divide -- fix the
    //      solution, never the oracle. ----
    {
        cpf::Rng rng(0xC0FFEEuLL);
        int cases = 0, diffs = 0;
        int first_bad_n = -1, first_bad_m = -1, first_bad_step = -1;
        constexpr int kIters = 6000;

        for (int t = 0; t < kIters; ++t) {
            const int n = static_cast<int>(rng.in_range(1, 64));
            const int m = static_cast<int>(rng.in_range(1, 64));
            std::vector<i64> A = rand_poly(n, rng);
            std::vector<i64> pts = rand_points(m, rng);

            // sprinkle x=0 and repeats into some cases -- the corners Horner and
            // the tree are most likely to disagree on.
            if (t % 3 == 0 && !pts.empty()) pts[0] = 0;
            if (t % 5 == 0 && pts.size() >= 2) pts[1] = pts[0];

            const std::vector<i64> got = multipoint_eval(A, pts);
            const std::vector<i64> want = horner_multipoint(A, pts);

            ++cases;
            if (got != want) {
                if (diffs == 0) {
                    first_bad_n = n;
                    first_bad_m = m;
                    first_bad_step = t;
                }
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs)
            std::printf("  first bad case: n=%d m=%d (generator step %d)\n",
                        first_bad_n, first_bad_m, first_bad_step);
        std::printf("differential: %d instances (n,m in 1..64), %d diffs\n",
                    cases, diffs);
    }

    // ---- a few larger sizes vs Horner, incl. 2^k and 2^k +/- 1, to exercise
    //      multi-level trees and the power-of-two padding in convolve. ----
    {
        cpf::Rng rng(0xBADC0DEuLL);
        int diffs = 0;
        for (int sz : {100, 255, 256, 257, 511, 512, 1000, 1024, 2000}) {
            std::vector<i64> A = rand_poly(sz, rng);
            std::vector<i64> pts = rand_points(sz, rng);
            if (multipoint_eval(A, pts) != horner_multipoint(A, pts)) ++diffs;
            // and the asymmetric shapes at scale.
            std::vector<i64> A2 = rand_poly(sz, rng);
            std::vector<i64> pts2 = rand_points(sz / 4 + 1, rng);
            if (multipoint_eval(A2, pts2) != horner_multipoint(A2, pts2)) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("larger sizes 100..2000 (incl. 2^k and 2^k +/- 1): %d diffs\n",
                    diffs);
    }

    return cpf::report();
}
