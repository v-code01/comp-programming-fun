#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0037::i64;
using p0037::u64;
using p0037::kMod;

namespace {

// build the first m terms of a sequence from a known recurrence + seeds. this is
// the generator's ground truth -- the exact bytes the solver will see, produced
// by the same naive definition the oracle trusts.
std::vector<i64> emit(const std::vector<i64>& rec, const std::vector<i64>& init,
                      int m) {
    std::vector<i64> a(init);
    a.resize(static_cast<std::size_t>(m));
    const int k = static_cast<int>(rec.size());
    for (int i = k; i < m; ++i) {
        __int128 acc = 0;
        for (int j = 0; j < k; ++j)
            acc += static_cast<__int128>(rec[j]) * a[i - 1 - j];
        a[static_cast<std::size_t>(i)] = static_cast<i64>(acc % kMod);
    }
    return a;
}

}  // namespace

int main() {
    // ---- hand examples: sequences whose n-th term everyone already knows. ----

    // fibonacci 0,1,1,2,3,5,8,13. order-2 law a_i = a_{i-1}+a_{i-2}.
    {
        std::vector<i64> fib{0, 1, 1, 2, 3, 5, 8, 13};
        CPF_EQ(p0037::solve(10, fib), 55);     // F(10) = 55.
        CPF_EQ(p0037::solve(20, fib), 6765);   // F(20) = 6765.
        CPF_EQ(p0037::solve(7, fib), 13);      // N < m -- read straight off.
        // BM must recover exactly (c_1,c_2) = (1,1).
        std::vector<i64> rec = p0037::berlekamp_massey(fib);
        CPF_EQ(static_cast<int>(rec.size()), 2);
        if (rec.size() == 2) {
            CPF_EQ(rec[0], 1);
            CPF_EQ(rec[1], 1);
        }
    }

    // geometric 1,2,4,8. order-1 law a_i = 2 a_{i-1}, so a_N = 2^N mod p.
    {
        std::vector<i64> geo{1, 2, 4, 8};
        for (u64 N : {u64(5), u64(10), u64(30), u64(100), u64(1000000000000000000ULL)}) {
            i64 want = p0037::pw(2, static_cast<i64>(N));
            CPF_EQ(p0037::solve(N, geo), want);
        }
        std::vector<i64> rec = p0037::berlekamp_massey(geo);
        CPF_EQ(static_cast<int>(rec.size()), 1);
        if (rec.size() == 1) CPF_EQ(rec[0], 2);
    }

    // constant 7,7,7,7. order-1 law a_i = a_{i-1}, so a_N = 7 for every N.
    {
        std::vector<i64> con{7, 7, 7, 7};
        CPF_EQ(p0037::solve(0, con), 7);
        CPF_EQ(p0037::solve(3, con), 7);
        CPF_EQ(p0037::solve(1000, con), 7);
        CPF_EQ(p0037::solve(1000000000000000000ULL, con), 7);
    }

    // ---- edges the code has to survive. ----

    // all-zero sequence: BM returns an empty recurrence, answer is 0 everywhere.
    {
        std::vector<i64> z(6, 0);
        CPF_EQ(static_cast<int>(p0037::berlekamp_massey(z).size()), 0);
        CPF_EQ(p0037::solve(0, z), 0);
        CPF_EQ(p0037::solve(5, z), 0);                       // N < m branch.
        CPF_EQ(p0037::solve(1000000000000000000ULL, z), 0);  // k==0 branch.
    }

    // order-1 with a nontrivial multiplier and seed: a = 3, 15, 75, 375 (x5).
    {
        std::vector<i64> g{3, 15, 75, 375};
        std::vector<i64> rec = p0037::berlekamp_massey(g);
        CPF_EQ(static_cast<int>(rec.size()), 1);
        if (rec.size() == 1) CPF_EQ(rec[0], 5);
        // a_N = 3 * 5^N.
        CPF_EQ(p0037::solve(40, g), p0037::mul(3, p0037::pw(5, 40)));
        CPF_EQ(p0037::solve(1000000000000000000ULL, g),
               p0037::mul(3, p0037::pw(5, 1000000000000000000LL)));
    }

    // N == 0 across a random-ish order-3 sequence: must return a_0 untouched.
    {
        std::vector<i64> rec{5, 7, 3}, init{2, 9, 4};
        std::vector<i64> a = emit(rec, init, 12);
        CPF_EQ(p0037::solve(0, a), 2);
        CPF_EQ(p0037::solve(1, a), 9);
        // and it agrees with both oracles at N=0.
        CPF_EQ(p0037::solve(0, a), p0037::companion_matrix_power(rec, init, 0));
        CPF_EQ(p0037::solve(0, a), p0037::naive_iterate(rec, init, 0));
    }

    // ---- randomized differential vs BOTH oracles, thousands of cases. ----
    //
    // each case: draw an order-k recurrence (k in 1..8) and k seeds, emit the
    // first m = 2k+3 terms, then hammer the solver at small / medium / 1e18 N.
    // small+medium N also check naive_iterate; every N checks companion power.
    // the generated recurrence is the independent ground truth BM must rebuild.
    {
        cpf::Rng rng(0xC0FFEEuLL);
        int cases = 0, diffs = 0;
        u64 first_bad_seed_case = 0;
        constexpr int kIters = 4000;

        for (int t = 0; t < kIters; ++t) {
            const int k = static_cast<int>(rng.in_range(1, 8));
            std::vector<i64> rec(k), init(k);
            // draw coefficients; force c_k != 0 so the order is genuinely k and
            // the companion matrix is nonsingular (a clean order-k generator).
            for (int j = 0; j < k; ++j) rec[j] = rng.in_range(0, kMod - 1);
            if (rec[k - 1] == 0) rec[k - 1] = 1;
            for (int j = 0; j < k; ++j) init[j] = rng.in_range(0, kMod - 1);

            const int m = 2 * k + 3;
            std::vector<i64> a = emit(rec, init, m);

            // assemble the N values to probe.
            std::vector<u64> Ns;
            Ns.push_back(static_cast<u64>(rng.in_range(0, m - 1)));   // small, in-prefix.
            Ns.push_back(static_cast<u64>(m));                        // first leap.
            Ns.push_back(static_cast<u64>(rng.in_range(m, 100000)));  // medium.
            Ns.push_back(1000000000000000000ULL);                    // the hard one.
            Ns.push_back(static_cast<u64>(rng.in_range(1000000000000LL,
                                                       999999999999999999LL)));

            for (u64 N : Ns) {
                const i64 got = p0037::solve(N, a);
                const i64 want_c = p0037::companion_matrix_power(rec, init, N);
                ++cases;
                bool ok = (got == want_c);
                if (N <= 100000) {  // reachable by honest iteration too.
                    const i64 want_n = p0037::naive_iterate(rec, init, N);
                    ok = ok && (got == want_n);
                }
                if (!ok) {
                    if (diffs == 0) first_bad_seed_case = static_cast<u64>(t);
                    ++diffs;
                }
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs)
            std::printf("  first diff at generator case %llu\n",
                        static_cast<unsigned long long>(first_bad_seed_case));
        std::printf("differential: %d checks over %d generated recurrences, "
                    "%d diffs\n",
                    cases, kIters, diffs);
    }

    return cpf::report();
}
