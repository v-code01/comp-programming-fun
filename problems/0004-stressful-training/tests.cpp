#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0004::min_charger_power;
using V = std::vector<std::int64_t>;

int main() {
    // the four statement examples -- ground truth.
    CPF_EQ(min_charger_power(V{3, 2}, V{4, 2}, 4), 5);   // needs power 5.
    CPF_EQ(min_charger_power(V{4}, V{2}, 5), 1);         // one student, k=5.
    CPF_EQ(min_charger_power(V{4}, V{2}, 6), 2);         // same but a minute longer.
    CPF_EQ(min_charger_power(V{2, 10}, V{3, 15}, 2), -1);  // two die at once, one step.

    // one student who never dips -- a huge head start, tiny drain. x = 0 works.
    CPF_EQ(min_charger_power(V{1000000000000LL}, V{1}, 2), 0);

    // nobody drains negative on their own, so the charger sits idle. x = 0.
    CPF_EQ(min_charger_power(V{100, 100}, V{1, 1}, 2), 0);

    // two students both hit zero-runway at the very first beginning, but k = 2
    // gives a single charge -- one lives, one dies, at any power. -1.
    CPF_EQ(min_charger_power(V{1, 1}, V{2, 2}, 2), -1);

    // overflow guard: b near 1e12 so thr = b*(k-1) and the num staircase run
    // large. answers hand-derived, so a silent int64 wrap would show as a diff.
    //   n=1, a=1, b=1e12, k=2: one step, need x >= 1e12-1.
    CPF_EQ(min_charger_power(V{1}, V{1000000000000LL}, 2), 999999999999LL);
    //   n=1, a=5, b=1e12, k=3: two steps, the beginning-3 demand binds at
    //   ceil((2e12-5)/x) <= 2, i.e. x >= 1e12-2.
    CPF_EQ(min_charger_power(V{5}, V{1000000000000LL}, 3), 999999999998LL);

    // single-student sanity across a few k: min x rises as the contest runs long.
    CPF_EQ(min_charger_power(V{4}, V{2}, 2), 0);  // a=4 >= b*(k-1)=2, no charge.
    CPF_EQ(min_charger_power(V{4}, V{2}, 3), 0);  // a=4 >= b*(k-1)=4, still fine.
    CPF_EQ(min_charger_power(V{4}, V{2}, 4), 1);  // now a < b*(k-1)=6, needs 1.

    // differential vs the greedy oracle on thousands of small boards, infeasible
    // ones included -- small a, big b, short k make plenty of -1s. one struct in,
    // one int64 out; a disagreement prints the seed that repro's it.
    struct Board {
        V a, b;
        std::int64_t k;
    };
    bool ok = cpf::differential(
        20000, 0xC0FFEE,
        [](cpf::Rng& rng) {
            std::int64_t n = rng.in_range(1, 6);
            std::int64_t k = rng.in_range(2, 8);
            V a(static_cast<std::size_t>(n)), b(static_cast<std::size_t>(n));
            for (auto& x : a) x = rng.in_range(1, 12);
            for (auto& x : b) x = rng.in_range(1, 12);
            return Board{a, b, k};
        },
        [](const Board& in) { return min_charger_power(in.a, in.b, in.k); },
        [](const Board& in) {
            return p0004::ReferenceSolver{}.solve(in.a, in.b, in.k);
        });
    CPF_CHECK(ok);

    // widen the value range so ceils and deadlines spread out more.
    bool ok_wide = cpf::differential(
        20000, 0x5EED,
        [](cpf::Rng& rng) {
            std::int64_t n = rng.in_range(1, 4);
            std::int64_t k = rng.in_range(2, 12);
            V a(static_cast<std::size_t>(n)), b(static_cast<std::size_t>(n));
            for (auto& x : a) x = rng.in_range(1, 40);
            for (auto& x : b) x = rng.in_range(1, 40);
            return Board{a, b, k};
        },
        [](const Board& in) { return min_charger_power(in.a, in.b, in.k); },
        [](const Board& in) {
            return p0004::ReferenceSolver{}.solve(in.a, in.b, in.k);
        });
    CPF_CHECK(ok_wide);

    std::printf("differential: 40000 random boards vs greedy oracle\n");
    return cpf::report();
}
