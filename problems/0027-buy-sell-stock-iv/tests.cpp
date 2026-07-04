#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0027::brute;
using p0027::i64;
using p0027::reference;
using p0027::solve;
using Vec = std::vector<int>;

namespace {

// the unlimited-trades profit -- sum of every positive daily move. the aliens
// answer must hit this the instant K clears the trade count it takes.
i64 greedy_unlimited(const Vec& p) {
    i64 s = 0;
    for (std::size_t i = 1; i < p.size(); ++i)
        if (p[i] > p[i - 1]) s += p[i] - p[i - 1];
    return s;
}

}  // namespace

int main() {
    // --- hand examples from the statement. ---
    // [3,2,6,5,0,3], K=2 -> buy2 sell6 (+4), buy0 sell3 (+3) = 7.
    CPF_EQ(solve(2, Vec{3, 2, 6, 5, 0, 3}), i64{7});
    // same prices, one trade -> the single best move, buy2 sell6 = 4.
    CPF_EQ(solve(1, Vec{3, 2, 6, 5, 0, 3}), i64{4});
    // K = 0 -> no trades allowed, profit 0.
    CPF_EQ(solve(0, Vec{3, 2, 6, 5, 0, 3}), i64{0});
    // classic [7,1,5,3,6,4], K=2 -> buy1 sell5 (+4), buy3 sell6 (+3) = 7.
    CPF_EQ(solve(2, Vec{7, 1, 5, 3, 6, 4}), i64{7});
    // [1,2,3,4,5] one long climb, K=2 -> one trade grabs it all = 4.
    CPF_EQ(solve(2, Vec{1, 2, 3, 4, 5}), i64{4});

    // --- K >= n/2 collapses to unlimited: greedy sum of positive moves. ---
    {
        Vec p{1, 7, 2, 8, 3, 9, 4, 10};
        i64 want = greedy_unlimited(p);  // 6 + 6 + 6 + 6 = 24.
        CPF_EQ(solve(4, p), want);       // K = 4 = n/2 -- exactly enough.
        CPF_EQ(solve(100, p), want);     // K way past n -- still unlimited.
        CPF_EQ(want, i64{24});
    }

    // --- edges. ---
    CPF_EQ(solve(5, Vec{4}), i64{0});          // n = 1 -- no later day to sell.
    CPF_EQ(solve(0, Vec{4}), i64{0});          // n = 1, K = 0.
    CPF_EQ(solve(3, Vec{1, 2, 3, 4, 5}), i64{4});      // monotone up.
    CPF_EQ(solve(3, Vec{5, 4, 3, 2, 1}), i64{0});      // monotone down -- sit out.
    CPF_EQ(solve(2, Vec{5, 5, 5, 5}), i64{0});         // flat -- no move to catch.
    CPF_EQ(solve(9, Vec{7, 7, 7, 7, 7}), i64{0});      // all equal, big K.
    // zero moves are free to hold across: 1 2 2 3 is ONE trade worth 2, not two.
    CPF_EQ(solve(1, Vec{1, 2, 2, 3}), i64{2});
    // a loss move splits it: need two trades to skip the dip, one trade caps at 1.
    CPF_EQ(solve(1, Vec{1, 2, 2, 1, 2}), i64{1});
    CPF_EQ(solve(2, Vec{1, 2, 2, 1, 2}), i64{2});

    // --- int64 headroom. big prices, many trades -- profit blows past 2^31. ---
    {
        Vec p;
        p.reserve(2000);
        for (int i = 0; i < 1000; ++i) { p.push_back(1); p.push_back(1000000000); }
        i64 want = static_cast<i64>(1000) * (1000000000 - 1);  // 999,999,999,000.
        CPF_EQ(solve(1000, p), want);
        CPF_EQ(want, i64{999999999000LL});
    }

    // --- the oracle itself, cross-checked against the exponential brute. ---
    // tiny n so brute's 3^n stays cheap. if the O(n*K) dp is wrong, this catches
    // it before it's ever used as a reference.
    {
        int cases = 0;
        for (std::uint64_t seed = 0; seed < 3000; ++seed) {
            cpf::Rng r(0xA11CE5 + seed);
            int n = static_cast<int>(r.in_range(1, 8));
            Vec p(static_cast<std::size_t>(n));
            for (auto& x : p) x = static_cast<int>(r.in_range(1, 6));
            for (int K = 0; K <= n; ++K) {
                ++cases;
                CPF_EQ(reference(K, p), brute(K, p));
            }
        }
        std::printf("oracle vs brute: %d (n,K,prices) cases\n", cases);
    }

    // --- the main event: aliens vs the O(n*K) oracle, ALL K from 0..n, on
    // thousands of small arrays. the shapes that stress trade counting -- flat,
    // monotone both ways, zigzag, tight ranges that force plateaus -- all fall
    // out of the random draws, and each array is swept over every K. ---
    {
        int arrays = 0, cases = 0, diffs = 0;
        int fseed = -1, fK = -1;
        // four regimes: wide range, tight range (many ties/plateaus), tiny
        // range (lots of flats), and longer arrays.
        const int lo_hi[4][3] = {
            {1, 10, 8},   // n<=8, prices 1..10.
            {1, 3, 10},   // n<=10, prices 1..3  -- plateaus everywhere.
            {1, 2, 12},   // n<=12, prices 1..2  -- flats and single moves.
            {1, 1000, 9}, // n<=9, prices 1..1000 -- wide, distinct slopes.
        };
        for (int regime = 0; regime < 4; ++regime) {
            const int plo = lo_hi[regime][0];
            const int phi = lo_hi[regime][1];
            const int nmax = lo_hi[regime][2];
            for (std::uint64_t seed = 0; seed < 4000; ++seed) {
                cpf::Rng r(0xBADC0FFEE0DDULL + regime * 1000000ULL + seed);
                int n = static_cast<int>(r.in_range(1, nmax));
                Vec p(static_cast<std::size_t>(n));
                for (auto& x : p) x = static_cast<int>(r.in_range(plo, phi));
                ++arrays;
                for (int K = 0; K <= n; ++K) {
                    ++cases;
                    i64 got = solve(K, p);
                    i64 want = reference(K, p);
                    if (got != want) {
                        if (diffs == 0) { fseed = static_cast<int>(seed); fK = K; }
                        ++diffs;
                    }
                }
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) std::printf("  first diff: regime seed %d, K %d\n", fseed, fK);
        std::printf("aliens vs oracle, all K: %d arrays, %d cases, %d diffs\n",
                    arrays, cases, diffs);
    }

    // --- explicit monotone / flat / zigzag sweeps over every K, hammered. ---
    {
        auto sweep = [](const Vec& p) {
            for (int K = 0; K <= static_cast<int>(p.size()); ++K)
                CPF_EQ(solve(K, p), reference(K, p));
        };
        sweep(Vec{1, 2, 3, 4, 5, 6, 7, 8});          // strictly up.
        sweep(Vec{8, 7, 6, 5, 4, 3, 2, 1});          // strictly down.
        sweep(Vec{5, 5, 5, 5, 5, 5});                // flat.
        sweep(Vec{1, 9, 1, 9, 1, 9, 1, 9});          // hard zigzag.
        sweep(Vec{3, 3, 5, 5, 2, 2, 8, 8, 1, 1});    // zigzag with plateaus.
        sweep(Vec{1, 2, 2, 3, 3, 2, 2, 4, 4});       // zeros between gains.
        std::printf("shape sweeps: passed\n");
    }

    return cpf::report();
}
