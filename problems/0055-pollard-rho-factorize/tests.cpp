#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0055::u64;
using u128 = unsigned __int128;

namespace {

// a prime in [lo, hi], drawn by rejection off the fast primality test. used to
// build large semiprimes the trial oracle can't reach.
u64 rand_prime(cpf::Rng& rng, u64 lo, u64 hi) {
    while (true) {
        u64 x = static_cast<u64>(rng.in_range(static_cast<std::int64_t>(lo),
                                              static_cast<std::int64_t>(hi))) |
                1ull;
        if (p0055::is_prime(x)) return x;
    }
}

// second, independent primality gate for the 1e18 self-check: > 1 and no divisor
// up to min(sqrt f, 1e6). below 1e12 that's a complete test; above it, a partial
// one that still catches any small factor rho might have left behind.
bool no_small_divisor(u64 f) {
    if (f < 2) return false;
    const u64 lim = 1000000;
    for (u64 d = 2; d <= lim && d <= f / d; ++d) {
        if (f % d == 0) return false;
    }
    return true;
}

}  // namespace

int main() {
    using namespace p0055;

    // ---- hand examples from the statement ----
    CPF_EQ(factorize(12), (std::vector<u64>{2, 2, 3}));  // verified by hand.
    CPF_EQ(factorize(999999999989ull),
           (std::vector<u64>{999999999989ull}));  // a prime -> itself.
    {
        std::vector<u64> fifty_nine_twos(59, 2);
        CPF_EQ(factorize(1ull << 59), fifty_nine_twos);  // 2^59 -> fifty-nine 2s.
    }
    CPF_EQ(factorize(1).empty(), true);  // 1 -> empty list.

    // ---- edges ----
    CPF_EQ(factorize(1).empty(), true);
    CPF_EQ(factorize(2), (std::vector<u64>{2}));
    CPF_EQ(factorize(1000003ull * 1000003ull),
           (std::vector<u64>{1000003ull, 1000003ull}));  // p^2, a rho stall risk.
    {
        // product of two ~1e9 primes -- a big semiprime near 1e18.
        u64 p = 999999937ull, q = 1000000007ull;  // both prime.
        std::vector<u64> want = {p, q};  // p < q, already sorted.
        CPF_EQ(factorize(p * q), want);
    }
    // largest prime below 1e18 range still round-trips to itself.
    CPF_EQ(is_prime(999999999999999989ull), true);
    CPF_EQ(factorize(999999999999999989ull),
           (std::vector<u64>{999999999999999989ull}));

    // ---- differential vs trial division, n <= 1e12 (thousands, mixed shapes) ----
    cpf::Rng rng(20250713);
    int diff_cases = 0, diffs = 0;
    u64 first_bad = 0;

    auto check_against_trial = [&](u64 n) {
        ++diff_cases;
        if (factorize(n) != trial_factor(n)) {
            if (diffs == 0) first_bad = n;
            ++diffs;
        }
    };

    // random uniform in [1, 1e12].
    for (int i = 0; i < 4000; ++i) {
        check_against_trial(
            static_cast<u64>(rng.in_range(1, 1000000000000ll)));
    }
    // structured: primes themselves (near 1e6, the rho worst-case scale).
    for (int i = 0; i < 1000; ++i) {
        check_against_trial(rand_prime(rng, 100000, 1000000));
    }
    // structured: semiprimes p*q with both factors near 1e6 -- one big-ish factor
    // to find, rho's hard case, product still under the oracle's 1e12 reach.
    for (int i = 0; i < 2000; ++i) {
        u64 p = rand_prime(rng, 100000, 1000000);
        u64 q = rand_prime(rng, 100000, 1000000);
        check_against_trial(p * q);
    }
    // structured: prime powers p^k.
    for (int i = 0; i < 1000; ++i) {
        u64 p = rand_prime(rng, 2, 40000);
        u64 v = p;
        while (v <= 1000000000000ull / p) v *= p;  // largest p^k <= 1e12.
        check_against_trial(v);
    }
    // structured: highly composite / smooth numbers -- many small factors.
    for (u64 base : {720720ull, 735134400ull, 999999999999ull, 1000000000000ull,
                     6469693230ull, 2ull * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23}) {
        check_against_trial(base);
    }
    // structured: 1.
    check_against_trial(1);

    CPF_CHECK(diffs == 0);
    if (diffs) {
        std::printf("  first differential mismatch at n = %llu\n",
                    static_cast<unsigned long long>(first_bad));
    }
    std::printf("differential vs trial (n<=1e12): %d cases, %d diffs\n",
                diff_cases, diffs);

    // ---- self-consistency on the FULL 1e18 range (trial can't reach here) ----
    // for each n: (a) factors multiply back to n exactly (u128), (b) every factor
    // clears an independent primality gate, (c) the list is non-decreasing.
    int self_cases = 0, self_bad = 0;
    u64 self_first_bad = 0;

    auto self_check = [&](u64 n) {
        ++self_cases;
        auto f = factorize(n);
        bool ok = true;
        u128 prod = 1;
        u64 prev = 0;
        for (u64 x : f) {
            prod *= x;
            if (x < 2 || !no_small_divisor(x)) ok = false;  // (b)
            if (x < prev) ok = false;                       // (c)
            prev = x;
        }
        if (prod != static_cast<u128>(n)) ok = false;  // (a)
        if (n == 1 && !f.empty()) ok = false;
        if (!ok) {
            if (self_bad == 0) self_first_bad = n;
            ++self_bad;
        }
    };

    // random uniform across the whole range.
    for (int i = 0; i < 4000; ++i) {
        self_check(static_cast<u64>(rng.in_range(1, 1000000000000000000ll)));
    }
    // structured large semiprimes: two ~1e9 primes, product ~1e18. this is the
    // regime rho is built for and trial division cannot touch.
    for (int i = 0; i < 2000; ++i) {
        u64 p = rand_prime(rng, 900000000, 1000000000);
        u64 q = rand_prime(rng, 900000000, 1000000000);
        self_check(p * q);
    }
    // structured: large prime powers and a lone large prime.
    for (int i = 0; i < 500; ++i) {
        self_check(rand_prime(rng, 1000000000000000000ull / 2,
                              1000000000000000000ull));
    }

    CPF_CHECK(self_bad == 0);
    if (self_bad) {
        std::printf("  first self-check failure at n = %llu\n",
                    static_cast<unsigned long long>(self_first_bad));
    }
    std::printf("self-check multiply-back + primality (n<=1e18): %d cases, "
                "%d failures\n",
                self_cases, self_bad);

    return cpf::report();
}
