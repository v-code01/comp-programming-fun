#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// independent 2^e mod p by fast exponentiation -- lets the mod-boundary test
// check solve() against a value the SOS path never touches.
std::int64_t pow2_mod(std::int64_t e) {
    std::int64_t base = 2, r = 1;
    while (e > 0) {
        if (e & 1) r = r * base % p0007::kMod;
        base = base * base % p0007::kMod;
        e >>= 1;
    }
    return r;
}

// one differential case: n random elements over a random small bit-width, kept
// small so AND == 0 actually happens and the 2^n brute stays cheap.
struct Case {
    std::vector<std::uint32_t> a;
    int bits;
};

}  // namespace

int main() {
    using namespace p0007;

    // the three statement examples.
    CPF_EQ(solve({2, 3, 3}, 20), 0);
    CPF_EQ(solve({0, 1, 2, 3}, 20), 10);
    CPF_EQ(solve({5, 2, 0, 5, 2, 1}, 20), 53);

    // edges the formula has to get right.
    CPF_EQ(solve({0}, 20), 1);          // single zero -- the group {0} works.
    CPF_EQ(solve({7}, 20), 0);          // single non-zero -- its AND is itself.
    CPF_EQ(solve({0, 0, 0}, 20), 7);    // all zero -- every non-empty group, 2^3-1.
    CPF_EQ(solve({1, 1, 1}, 20), 0);    // shared bit never clears -- no group.
    CPF_EQ(solve({1, 2}, 20), 1);       // disjoint bits -- only {1,2} ANDs to 0.

    // bits = 0 corner: the only value is 0, one mask, answer is 2^n - 1.
    CPF_EQ(solve({0, 0, 0, 0}, 0), 15);

    // mod-boundary: n zeros gives exactly 2^n - 1, big enough to wrap the modulus.
    // check against an independent modpow, not the brute -- brute can't reach n here.
    for (std::int64_t n : {30, 100000, 1000000}) {
        std::vector<std::uint32_t> zeros(static_cast<std::size_t>(n), 0);
        std::int64_t want = (pow2_mod(n) - 1 + kMod) % kMod;
        CPF_EQ(solve(zeros, 20), want);
    }

    // randomized differential vs the 2^n brute. thousands of seeds, each a fresh
    // (n, bit-width, values) draw -- narrow bit-widths force real AND == 0 hits.
    bool ok = cpf::differential(
        20000, 0xC0FFEE,
        [](cpf::Rng& rng) {
            Case c;
            c.bits = static_cast<int>(rng.in_range(1, 6));
            int n = static_cast<int>(rng.in_range(1, 16));
            std::uint32_t hi = (1u << c.bits) - 1;
            c.a.resize(static_cast<std::size_t>(n));
            for (auto& x : c.a)
                x = static_cast<std::uint32_t>(rng.in_range(0, hi));
            return c;
        },
        [](const Case& c) { return solve(c.a, c.bits); },
        [](const Case& c) { return brute(c.a); });
    CPF_CHECK(ok);
    std::printf("differential: 20000 random cases vs 2^n brute\n");

    return cpf::report();
}
