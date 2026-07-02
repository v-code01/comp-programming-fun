#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0011::min_cost;
using V = std::vector<std::int64_t>;

// build a valid board straight from the guarantees: a strictly increasing with
// a_1 = 1, b strictly decreasing with b_n = 0, every b_i>0 before the last. steps
// are drawn in [1, step] so both stay strictly monotone and inside their bounds.
struct Board {
    V a, b;
};

static Board make_board(cpf::Rng& rng, std::int64_t n, std::int64_t step) {
    const std::size_t un = static_cast<std::size_t>(n);
    V a(un), b(un);
    a[0] = 1;
    for (std::size_t i = 1; i < un; ++i) a[i] = a[i - 1] + rng.in_range(1, step);
    b[un - 1] = 0;  // b_n = 0 by guarantee.
    for (std::size_t i = un - 1; i-- > 0;) b[i] = b[i + 1] + rng.in_range(1, step);
    return {a, b};
}

int main() {
    // the two statement examples -- ground truth.
    CPF_EQ(min_cost(V{1, 2, 3, 4, 5}, V{5, 4, 3, 2, 0}), 25);
    CPF_EQ(min_cost(V{1, 2, 3, 10, 20, 30}, V{6, 5, 4, 3, 2, 0}), 138);

    // n = 1: a single height-1 tree falls off the initial charge. free.
    CPF_EQ(min_cost(V{1}, V{0}), 0);

    // n = 2: tree 1 free, then tree 2 needs a_2 chops each costing b_1.
    //   dp[2] = a_2 * b_1 = 2 * 7 = 14.
    CPF_EQ(min_cost(V{1, 2}, V{7, 0}), 14);
    //   a taller second tree pays the same b_1 per chop: 100 * 5 = 500.
    CPF_EQ(min_cost(V{1, 100}, V{5, 0}), 500);

    // the base case where the hull only ever holds line j=0 until the end -- three
    // trees, and leaning on tree 1 straight through beats any relay.
    //   dp[2] = 2*9 = 18; dp[3] = min(3*9, 18 + 3*4) = min(27, 30) = 27.
    CPF_EQ(min_cost(V{1, 2, 3}, V{9, 4, 0}), 27);
    //   here the relay wins: dp[2] = 2*10 = 20; dp[3] = min(50*10, 20 + 50*1)
    //   = min(500, 70) = 70 -- felling tree 2 first then chopping is far cheaper.
    CPF_EQ(min_cost(V{1, 2, 50}, V{10, 1, 0}), 70);

    // overflow guard: values pinned near 1e9 so a_i*b_j lands at ~1e18 and the
    // hull's cross product at ~1e27. answers hand-derived, so an int64 wrap in
    // eval() or a non-__int128 bad() would show up as a mismatch here.
    //   n=2: dp[2] = a_2 * b_1 = 1000000000 * 1000000000 = 1e18.
    CPF_EQ(min_cost(V{1, 1000000000LL}, V{1000000000LL, 0}),
           1000000000000000000LL);
    //   n=3, all near 1e9. dp[2] = a_2*b_1 = 999999999 * 1000000000
    //   = 999999999000000000. dp[3] = min(a_3*b_1, dp[2] + a_3*b_2):
    //     a_3*b_1 = 1000000000 * 1000000000        = 1000000000000000000
    //     dp[2] + a_3*b_2 = 999999999000000000
    //                       + 1000000000 * 999999999
    //                     = 999999999000000000 + 999999999000000000
    //                     = 1999999998000000000  -> first term wins.
    CPF_EQ(min_cost(V{1, 999999999LL, 1000000000LL},
                    V{1000000000LL, 999999999LL, 0}),
           1000000000000000000LL);

    // the CHT solution must match the O(n^2) oracle everywhere. thousands of
    // random valid boards, small so the oracle stays fast; a disagreement prints
    // the seed that replays it.
    bool ok_small = cpf::differential(
        20000, 0x319C,
        [](cpf::Rng& rng) {
            return make_board(rng, rng.in_range(1, 8), 6);
        },
        [](const Board& in) { return min_cost(in.a, in.b); },
        [](const Board& in) {
            return p0011::reference_min_cost(in.a, in.b);
        });
    CPF_CHECK(ok_small);

    // wider boards and bigger steps so the hull grows several lines deep and the
    // pointer actually walks -- exercises pops that retire the resting line.
    bool ok_wide = cpf::differential(
        20000, 0xD119A,
        [](cpf::Rng& rng) {
            return make_board(rng, rng.in_range(2, 40), 50);
        },
        [](const Board& in) { return min_cost(in.a, in.b); },
        [](const Board& in) {
            return p0011::reference_min_cost(in.a, in.b);
        });
    CPF_CHECK(ok_wide);

    // large-value boards: steps near 1e8 push a and b across their full range, so
    // the differential also stresses the int64/__int128 arithmetic, not just the
    // hull logic.
    bool ok_big = cpf::differential(
        20000, 0xB16A,
        [](cpf::Rng& rng) {
            return make_board(rng, rng.in_range(2, 25), 100000000LL);
        },
        [](const Board& in) { return min_cost(in.a, in.b); },
        [](const Board& in) {
            return p0011::reference_min_cost(in.a, in.b);
        });
    CPF_CHECK(ok_big);

    std::printf("differential: 60000 random valid boards vs O(n^2) oracle\n");
    return cpf::report();
}
