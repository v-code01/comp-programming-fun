#include <algorithm>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0008::reference;
using p0008::solve;

// one random small case: n tiny so the O(n^2 k) oracle is instant, a small
// alphabet so distinct-counts actually vary and the window add/remove gets
// stressed, k anywhere in [1, n] so every piece count shows up.
struct Case {
    std::vector<int> a;
    int k;
};

static Case gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 12));
    int cap = std::min(n, 4);  // types stay in [1, n] as the constraints demand.
    int alpha = static_cast<int>(rng.in_range(1, cap));
    std::vector<int> a(static_cast<std::size_t>(n));
    for (auto& x : a) x = static_cast<int>(rng.in_range(1, alpha));
    int k = static_cast<int>(rng.in_range(1, n));
    return {a, k};
}

int main() {
    // the three statement examples.
    CPF_EQ(solve({1, 2, 2, 1}, 1), 2);
    CPF_EQ(solve({1, 3, 3, 1, 4, 4, 4}, 2), 5);
    CPF_EQ(solve({7, 7, 8, 7, 7, 8, 1, 7}, 3), 6);

    // edges the recurrence and the window have to get right:
    // k == 1 -- one piece is the whole row, worth its distinct count.
    CPF_EQ(solve({1, 2, 2, 1}, 1), 2);
    CPF_EQ(solve({4, 4, 4, 4, 4}, 1), 1);
    // k == n -- every cake its own piece, each worth 1, sum is n.
    CPF_EQ(solve({1, 2, 2, 1}, 4), 4);
    CPF_EQ(solve({2, 2, 2}, 3), 3);
    // all one type -- k pieces, each distinct count 1, total k.
    CPF_EQ(solve({5, 5, 5, 5, 5}, 3), 3);
    CPF_EQ(solve({5, 5, 5, 5, 5}, 1), 1);
    // all distinct types -- no repeats, so any split sums to n.
    CPF_EQ(solve({1, 2, 3, 4, 5}, 2), 5);
    CPF_EQ(solve({1, 2, 3, 4, 5}, 3), 5);
    // single cake, single piece (n == 1 forces type 1).
    CPF_EQ(solve({1}, 1), 1);

    // the whole risk is the window add/remove staying in lockstep with the cut
    // it claims to cover -- so hammer it. thousands of random small cases, fast
    // solver against the brute oracle, any diff is a seed you can replay.
    bool ok = cpf::differential(
        30000, 0xBA'CE'11, gen, [](const Case& c) { return solve(c.a, c.k); },
        [](const Case& c) { return reference(c.a, c.k); });
    CPF_CHECK(ok);
    std::printf("differential: 30000 random cases vs O(n^2 k) oracle\n");

    return cpf::report();
}
