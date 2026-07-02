#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0003::brute;
using p0003::solve;

using Ranges = std::vector<std::pair<int, int>>;

// one random small case: n and q small enough that the O(q^2*n) oracle is
// instant, big enough that cover 1/2/3+ all show up and the pair math gets
// stressed.
struct Case {
    int n;
    int q;
    Ranges ranges;
    bool operator==(const Case&) const = default;
};

static Case gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 12));
    int q = static_cast<int>(rng.in_range(2, 8));
    Ranges r(static_cast<std::size_t>(q));
    for (auto& p : r) {
        int l = static_cast<int>(rng.in_range(1, n));
        int hi = static_cast<int>(rng.in_range(1, n));
        if (l > hi) std::swap(l, hi);
        p = {l, hi};
    }
    return {n, q, r};
}

int main() {
    // the three statement examples.
    CPF_EQ(solve(7, 5, Ranges{{1, 4}, {4, 5}, {5, 6}, {6, 7}, {3, 5}}),
           std::int64_t{7});
    CPF_EQ(solve(4, 3, Ranges{{1, 1}, {2, 2}, {3, 4}}), std::int64_t{2});
    CPF_EQ(solve(4, 4, Ranges{{1, 1}, {2, 2}, {2, 3}, {3, 4}}),
           std::int64_t{3});

    // edges the pair math has to survive:
    // q == 2 -- fire both, nothing left, answer 0.
    CPF_EQ(solve(5, 2, Ranges{{1, 5}, {2, 3}}), std::int64_t{0});
    // every painter identical -- firing two of many still leaves the range lit.
    CPF_EQ(solve(6, 5, Ranges{{2, 5}, {2, 5}, {2, 5}, {2, 5}, {2, 5}}),
           std::int64_t{4});
    // disjoint singles -- firing two loses exactly those two sections.
    CPF_EQ(solve(4, 4, Ranges{{1, 1}, {2, 2}, {3, 3}, {4, 4}}),
           std::int64_t{2});
    // single section, everyone piled on it -- can't lose it with only two fires.
    CPF_EQ(solve(1, 4, Ranges{{1, 1}, {1, 1}, {1, 1}, {1, 1}}),
           std::int64_t{1});
    // single section, only two painters -- fire both, it goes dark.
    CPF_EQ(solve(1, 2, Ranges{{1, 1}, {1, 1}}), std::int64_t{0});
    // painter id 0 is the sole coverer somewhere -- S1 == 0 must still count.
    CPF_EQ(solve(3, 3, Ranges{{1, 1}, {2, 2}, {3, 3}}), std::int64_t{1});
    // a cover-2 pair involving id 0 -- exercises the quadratic with i == 0.
    CPF_EQ(solve(2, 3, Ranges{{1, 2}, {1, 2}, {2, 2}}), std::int64_t{2});

    // the whole risk lives in the pair recovery, so hammer it -- thousands of
    // random small cases, fast solver against the brute oracle, any diff is a
    // seed you can replay.
    bool ok = cpf::differential(
        20000, 0xC0FFEE, gen,
        [](const Case& c) { return solve(c.n, c.q, c.ranges); },
        [](const Case& c) { return brute(c.n, c.q, c.ranges); });
    CPF_CHECK(ok);
    std::printf("differential: 20000 random cases vs brute oracle\n");

    return cpf::report();
}
