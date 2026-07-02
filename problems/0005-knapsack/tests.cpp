#include <array>
#include <cstdint>
#include <cstdio>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0005::solve;
using Cnt = std::array<std::int64_t, 9>;

static std::int64_t total_weight(const Cnt& c) {
    std::int64_t s = 0;
    for (int i = 1; i <= 8; ++i) s += static_cast<std::int64_t>(i) * c[i];
    return s;
}

int main() {
    // --- the three statement examples, ground truth. ---
    CPF_EQ(solve(10, Cnt{0, 1, 2, 3, 4, 5, 6, 7, 8}), 10);
    CPF_EQ(solve(0, Cnt{0, 0, 0, 0, 0, 0, 0, 0, 0}), 0);
    CPF_EQ(solve(3, Cnt{0, 0, 4, 1, 0, 0, 9, 8, 3}), 3);

    // --- edges the branches have to nail. ---
    CPF_EQ(solve(0, Cnt{0, 5, 5, 5, 5, 5, 5, 5, 5}), 0);   // W=0, empty set.
    CPF_EQ(solve(100, Cnt{0, 0, 0, 0, 0, 0, 0, 0, 0}), 0);  // nothing to take.
    CPF_EQ(solve(1000, Cnt{0, 3, 0, 0, 0, 0, 0, 0, 0}),
           3);  // total 3 < W -- take all.
    CPF_EQ(solve(15, Cnt{0, 0, 0, 5, 0, 0, 0, 0, 0}),
           15);  // five 3s hit W = 15 exactly.
    CPF_EQ(solve(14, Cnt{0, 0, 0, 5, 0, 0, 0, 0, 0}),
           12);  // 3s can't make 14 -- best is 12.

    // single weight, huge count, small W: pure multiples of 8.
    CPF_EQ(solve(10, Cnt{0, 0, 0, 0, 0, 0, 0, 0, 1000000000000000LL}),
           8);  // one 8 fits under 10, two don't.
    // single weight, huge count, large W but total < W -- take the whole pile.
    CPF_EQ(solve(1000000000000000000LL,
                 Cnt{0, 0, 0, 0, 0, 0, 0, 0, 10000000000000000LL}),
           80000000000000000LL);  // 8 * 1e16.
    // single weight 8, huge count, large W with total > W: largest multiple of 8.
    // 5e16 is divisible by 8, so the answer is 5e16 exactly.
    CPF_EQ(solve(50000000000000000LL,
                 Cnt{0, 0, 0, 0, 0, 0, 0, 0, 10000000000000000LL}),
           50000000000000000LL);
    // and one off it: 5e16 + 3 is not a multiple of 8, so shave to 5e16.
    CPF_EQ(solve(50000000000000003LL,
                 Cnt{0, 0, 0, 0, 0, 0, 0, 0, 10000000000000000LL}),
           50000000000000000LL);

    // W = 1e18 with big counts -- total 3.6e17 < W, so the answer is the total.
    CPF_EQ(solve(1000000000000000000LL,
                 Cnt{0, 10000000000000000LL, 10000000000000000LL,
                     10000000000000000LL, 10000000000000000LL,
                     10000000000000000LL, 10000000000000000LL,
                     10000000000000000LL, 10000000000000000LL}),
           360000000000000000LL);

    // --- thousands of small random cases vs the copy-by-copy oracle. ---
    // small W, small counts, zeros welcome. spans total < W, total = W, total > W
    // by construction. any disagreement prints its seed and stops.
    {
        cpf::Rng rng(0xC0FFEE);
        int cases = 0, diffs = 0;
        std::int64_t bad_W = -1;
        Cnt bad_c{};
        for (int t = 0; t < 40000; ++t) {
            std::int64_t W = rng.in_range(0, 400);
            Cnt c{};
            for (int i = 1; i <= 8; ++i) c[i] = rng.in_range(0, 6);
            ++cases;
            std::int64_t got = solve(W, c);
            std::int64_t want = p0005::reference_small(W, c);
            if (got != want) {
                if (diffs == 0) {
                    bad_W = W;
                    bad_c = c;
                }
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first diff: W=%lld cnt=", static_cast<long long>(bad_W));
            for (int i = 1; i <= 8; ++i)
                std::printf("%lld ", static_cast<long long>(bad_c[i]));
            std::printf("\n");
        }
        std::printf("small differential: %d cases, %d diffs\n", cases, diffs);
    }

    // --- the large branch, vs the independent bitset oracle. ---
    // W past the 13440 cutoff forces solve() into residue+block mode, and the
    // small oracle can't reach here. counts are sized so the total clears W (else
    // it's the trivial take-all), and W stays small enough for the bitset dp.
    {
        cpf::Rng rng(0xBEEF01);
        int cases = 0, diffs = 0;
        std::int64_t bad_W = -1;
        Cnt bad_c{};
        for (int t = 0; t < 4000; ++t) {
            std::int64_t W = rng.in_range(13440, 300000);
            Cnt c{};
            // spread weight so the total tends to exceed W and residues vary.
            for (int i = 1; i <= 8; ++i) c[i] = rng.in_range(0, 12000);
            if (total_weight(c) <= W) continue;  // trivial branch, skip.
            ++cases;
            std::int64_t got = solve(W, c);
            std::int64_t want = p0005::reference_medium(W, c);
            if (got != want) {
                if (diffs == 0) {
                    bad_W = W;
                    bad_c = c;
                }
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first diff: W=%lld cnt=", static_cast<long long>(bad_W));
            for (int i = 1; i <= 8; ++i)
                std::printf("%lld ", static_cast<long long>(bad_c[i]));
            std::printf("\n");
        }
        std::printf("large differential: %d cases, %d diffs\n", cases, diffs);
    }

    // --- sparse counts at large W: one or two weight types, huge multiplicities.
    // still checked against the bitset oracle at medium W where residues get
    // chunky and the lightest-representative argument earns its keep.
    {
        cpf::Rng rng(0x5040);
        int cases = 0, diffs = 0;
        for (int t = 0; t < 2000; ++t) {
            std::int64_t W = rng.in_range(13440, 300000);
            Cnt c{};
            int a = static_cast<int>(rng.in_range(1, 8));
            int b = static_cast<int>(rng.in_range(1, 8));
            c[a] = rng.in_range(0, 100000);
            c[b] = rng.in_range(0, 100000);
            if (total_weight(c) <= W) continue;
            ++cases;
            if (solve(W, c) != p0005::reference_medium(W, c)) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("sparse differential: %d cases, %d diffs\n", cases, diffs);
    }

    // --- huge W invariants, where no oracle can follow. ---
    // two things must hold, both proven in solution.hpp:
    //   answer <= W always, and
    //   when the pile outweighs W, the answer sits within 7 of W -- a subset can't
    //   step over the (W-8, W] window since each item weighs at most 8.
    {
        cpf::Rng rng(0xA11CE);
        int cases = 0, bad = 0;
        for (int t = 0; t < 20000; ++t) {
            std::int64_t W = rng.in_range(0, 1000000000000000000LL);
            Cnt c{};
            for (int i = 1; i <= 8; ++i)
                c[i] = rng.in_range(0, 10000000000000000LL);
            std::int64_t ans = solve(W, c);
            std::int64_t sum = total_weight(c);
            ++cases;
            bool ok = ans >= 0 && ans <= W;
            if (sum <= W)
                ok = ok && (ans == sum);
            else
                ok = ok && (ans >= W - 7);  // the proven within-8 window.
            if (!ok) ++bad;
        }
        CPF_CHECK(bad == 0);
        std::printf("huge-W invariants: %d cases, %d violations\n", cases, bad);
    }

    return cpf::report();
}
