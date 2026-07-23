#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Vec = std::vector<std::uint32_t>;

// a constant array of value v over 2^k subsets -- the all-zero and all-one fixtures.
Vec constant_array(int k, std::uint32_t v) {
    return Vec(static_cast<std::size_t>(1) << k, v);
}

}  // namespace

int main() {
    using namespace p0056;

    // --- hand examples, worked out by hand in the header math ---

    // k=0: the only subset is the empty set, c[0] = a[0]*b[0].
    CPF_EQ(solve({7}, {9}, 0), Vec({63}));

    // k=1: c[0] = a0*b0, c[1] = a0*b1 + a1*b0.
    //   a=(2,3), b=(5,7): c0 = 10, c1 = 2*7 + 3*5 = 29.
    CPF_EQ(solve({2, 3}, {5, 7}, 1), Vec({10, 29}));

    // k=2, worked out term by term with a=(1,2,3,4), b=(5,6,7,8):
    //   c0 = 1*5 = 5
    //   c1 = a0*b1 + a1*b0 = 6 + 10 = 16
    //   c2 = a0*b2 + a2*b0 = 7 + 15 = 22
    //   c3 = a0*b3 + a1*b2 + a2*b1 + a3*b0 = 8 + 14 + 18 + 20 = 60
    CPF_EQ(solve({1, 2, 3, 4}, {5, 6, 7, 8}, 2), Vec({5, 16, 22, 60}));

    // --- edges ---

    // all-zero in, all-zero out at a few widths.
    for (int k = 0; k <= 6; ++k) {
        CPF_EQ(solve(constant_array(k, 0), constant_array(k, 0), k),
               constant_array(k, 0));
    }

    // all-one: c[S] = sum over T subset S of 1*1 = (number of subsets of S)
    //        = 2^popcount(S). verify the k=2 corner by hand -- c[3] must be 4.
    {
        Vec c = solve(constant_array(2, 1), constant_array(2, 1), 2);
        CPF_EQ(c, Vec({1, 2, 2, 4}));  // 2^0, 2^1, 2^1, 2^2.
        CPF_EQ(c[3], 4u);              // the by-hand check the brief asks for.
    }
    // all-one for wider k: entry S is 2^popcount(S), reduced mod p.
    for (int k = 0; k <= 10; ++k) {
        Vec c = solve(constant_array(k, 1), constant_array(k, 1), k);
        int diffs = 0;
        for (std::size_t S = 0; S < c.size(); ++S) {
            const int pc = __builtin_popcountll(static_cast<unsigned long long>(S));
            std::uint32_t want = 1;
            for (int t = 0; t < pc; ++t) want = want * 2u % kMod;
            if (c[S] != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
    }

    // near-modulus values -- exercises the u64 product % p path against the oracle.
    {
        const std::uint32_t m = kMod - 1;  // p-1, so (p-1)^2 == 1 mod p.
        Vec c = solve({m, m}, {m, m}, 1);
        // c0 = (p-1)^2 = 1; c1 = 2*(p-1)^2 = 2.
        CPF_EQ(c, Vec({1, 2}));
        CPF_EQ(c, reference({m, m}, {m, m}, 1));
    }

    // --- randomized differential vs the O(3^k) submask oracle ---
    //
    // two sweeps. a wide-and-shallow one for volume, and a narrow-and-deep one so
    // the rank bookkeeping is exercised at real width. any disagreement prints the
    // seed and the first (k, S) where they part.
    auto run_sweep = [](int count, std::int64_t k_lo, std::int64_t k_hi,
                        std::uint64_t base_seed, const char* label) {
        int cases = 0, diffs = 0;
        std::uint64_t first_bad_seed = 0;
        int first_k = -1;
        std::size_t first_S = 0;
        for (int t = 0; t < count; ++t) {
            const std::uint64_t seed = base_seed + static_cast<std::uint64_t>(t);
            cpf::Rng rng(seed);
            const int k = static_cast<int>(rng.in_range(k_lo, k_hi));
            const std::size_t size = static_cast<std::size_t>(1) << k;
            Vec a(size), b(size);
            for (auto& x : a)
                x = static_cast<std::uint32_t>(rng.in_range(0, kMod - 1));
            for (auto& x : b)
                x = static_cast<std::uint32_t>(rng.in_range(0, kMod - 1));

            const Vec got = solve(a, b, k);
            const Vec want = reference(a, b, k);
            ++cases;
            if (got != want) {
                if (diffs == 0) {
                    first_bad_seed = seed;
                    first_k = k;
                    for (std::size_t S = 0; S < size; ++S) {
                        if (got[S] != want[S]) {
                            first_S = S;
                            break;
                        }
                    }
                }
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  DIFF %s: seed %llu, k=%d, first S=%zu\n", label,
                        static_cast<unsigned long long>(first_bad_seed), first_k,
                        first_S);
        }
        std::printf("differential %-14s %6d cases, %d diffs (k in [%lld,%lld])\n",
                    label, cases, diffs, static_cast<long long>(k_lo),
                    static_cast<long long>(k_hi));
        return cases;
    };

    int total = 0;
    total += run_sweep(8000, 0, 8, 0xC0FFEEu, "wide");    // volume, small k.
    total += run_sweep(400, 9, 12, 0xBADF00Du, "deep");   // depth, up to 3^12.
    std::printf("differential total: %d random instances vs O(3^k) oracle\n",
                total);

    return cpf::report();
}
