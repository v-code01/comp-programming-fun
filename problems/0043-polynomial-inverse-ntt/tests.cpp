#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0043::i64;
using p0043::kMod;

namespace {

// draw a random polynomial of length n with a_0 != 0 -- the invertibility
// precondition. coefficients uniform in [0, p).
std::vector<i64> rand_poly(int n, cpf::Rng& rng) {
    std::vector<i64> a(static_cast<std::size_t>(n));
    for (auto& x : a) x = rng.in_range(0, kMod - 1);
    if (a[0] == 0) a[0] = rng.in_range(1, kMod - 1);  // force a_0 != 0.
    return a;
}

}  // namespace

int main() {
    using p0043::inverse;
    using p0043::direct_inverse;
    using p0043::verify_inverse;

    // ---- hand examples worked out by hand ----
    {
        // A = [1] is the constant 1. its inverse is 1.
        CPF_EQ(inverse({1}, 1), (std::vector<i64>{1}));

        // A = 1 - x = [1, p-1]. 1/(1-x) = 1 + x + x^2 + … -> all ones. verify by
        // hand: (1 - x)(1 + x + … + x^{k-1}) = 1 - x^k == 1 (mod x^k).
        {
            std::vector<i64> A{1, kMod - 1};
            std::vector<i64> want(5, 1);
            CPF_EQ(inverse(A, 5), want);
            CPF_EQ(inverse(A, 1), (std::vector<i64>{1}));
            std::vector<i64> want12(12, 1);
            CPF_EQ(inverse(A, 12), want12);
        }

        // A = [2] is the constant 2. inverse is 2^{-1} mod p.
        CPF_EQ(inverse({2}, 1), (std::vector<i64>{p0043::inv(2)}));

        // A = 1 + x = [1, 1]. 1/(1+x) = 1 - x + x^2 - x^3 + … -> alternating
        // 1, p-1, 1, p-1, … (−1 mod p is p−1).
        {
            std::vector<i64> A{1, 1};
            std::vector<i64> want{1, kMod - 1, 1, kMod - 1, 1, kMod - 1};
            CPF_EQ(inverse(A, 6), want);
        }

        // a small poly cross-checked against the direct oracle explicitly.
        {
            std::vector<i64> A{3, 1, 4, 1, 5, 9, 2, 6};
            const int n = static_cast<int>(A.size());
            CPF_EQ(inverse(A, n), direct_inverse(A, n));
            CPF_CHECK(verify_inverse(A, inverse(A, n), n));
        }
    }

    // ---- edges ----
    {
        // n = 1: the answer is a single coefficient, a_0^{-1}.
        {
            std::vector<i64> A{123456789};
            CPF_EQ(inverse(A, 1), (std::vector<i64>{p0043::inv(123456789)}));
        }

        // a_0 large (near p): inverse still exact, and A*B == 1.
        {
            std::vector<i64> A{kMod - 1, 7, 0, 3, 999};
            const int n = static_cast<int>(A.size());
            CPF_EQ(inverse(A, n), direct_inverse(A, n));
            CPF_CHECK(verify_inverse(A, inverse(A, n), n));
        }

        // a middle zero coefficient -- newton must not choke on it.
        {
            std::vector<i64> A{5, 0, 0, 7, 0, 11, 0, 0, 0, 13};
            const int n = static_cast<int>(A.size());
            CPF_EQ(inverse(A, n), direct_inverse(A, n));
            CPF_CHECK(verify_inverse(A, inverse(A, n), n));
        }

        // truncation shorter than deg A: B holds only the first n coefficients.
        {
            std::vector<i64> A{2, 3, 5, 7, 11, 13, 17, 19};
            for (int n = 1; n <= static_cast<int>(A.size()); ++n) {
                CPF_EQ(inverse(A, n), direct_inverse(A, n));
                CPF_CHECK(verify_inverse(A, inverse(A, n), n));
            }
        }

        // n larger than deg A: request more precision than A has terms.
        {
            std::vector<i64> A{4, 9};  // 4 + 9x, ask for 20 coefficients of 1/A.
            CPF_EQ(inverse(A, 20), direct_inverse(A, 20));
            CPF_CHECK(verify_inverse(A, inverse(A, 20), 20));
        }
    }

    // ---- randomized differential vs the direct O(n^2) oracle AND the A*B==1
    //      verifier, thousands of small polynomials. two independent checks per
    //      case. any mismatch is the NTT or a newton truncation length -- fix the
    //      solution, never the oracle. ----
    {
        cpf::Rng rng(0xC0FFEEuLL);
        int cases = 0, direct_diffs = 0, verify_fails = 0;
        int first_bad_n = -1;
        std::uint64_t first_bad_seed = 0;
        constexpr int kIters = 6000;

        for (int t = 0; t < kIters; ++t) {
            const int n = static_cast<int>(rng.in_range(1, 64));
            const std::vector<i64> A = rand_poly(n, rng);

            const std::vector<i64> got = inverse(A, n);
            const std::vector<i64> want = direct_inverse(A, n);

            ++cases;
            const bool matches_direct = (got == want);
            const bool passes_verify = verify_inverse(A, got, n);
            if (!matches_direct || !passes_verify) {
                if (direct_diffs == 0 && verify_fails == 0) {
                    first_bad_n = n;
                    first_bad_seed = static_cast<std::uint64_t>(t);
                }
                if (!matches_direct) ++direct_diffs;
                if (!passes_verify) ++verify_fails;
            }
        }
        CPF_CHECK(direct_diffs == 0);
        CPF_CHECK(verify_fails == 0);
        if (direct_diffs || verify_fails)
            std::printf("  first bad case: n=%d (generator step %llu)\n",
                        first_bad_n,
                        static_cast<unsigned long long>(first_bad_seed));
        std::printf("differential: %d polys (n in 1..64), %d direct diffs, "
                    "%d verify fails\n",
                    cases, direct_diffs, verify_fails);
    }

    // ---- a few larger sizes, verifier only (direct oracle is O(n^2), still
    //      cheap here, so run both). exercises multi-doubling newton and the
    //      power-of-two overshoot on non-power-of-two n. ----
    {
        cpf::Rng rng(0xBADC0DEuLL);
        int diffs = 0;
        for (int n : {100, 255, 256, 257, 511, 512, 1000, 1023, 2000}) {
            const std::vector<i64> A = rand_poly(n, rng);
            const std::vector<i64> got = inverse(A, n);
            if (got != direct_inverse(A, n)) ++diffs;
            if (!verify_inverse(A, got, n)) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("larger sizes 100..2000 (incl. 2^k and 2^k+/-1): %d diffs\n",
                    diffs);
    }

    return cpf::report();
}
