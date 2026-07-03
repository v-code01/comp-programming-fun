#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0020::Fact;
using p0020::reference;
using p0020::solve;
using Facts = std::vector<Fact>;

// what the differential feeds both solvers: array size + the fact list.
struct Input {
    int n;
    Facts facts;
};

int main() {
    // --- the three statement examples. each hand-confirmed with an independent
    //     brute before hard-coding, and cross-checked against the oracle here. ---
    CPF_EQ(solve(2, Facts{}), 2);                              // n=2, no facts.
    CPF_EQ(solve(3, Facts{{1, 1, 3, 2}}), 5);                 // all a_i >= 2.
    CPF_EQ(solve(3, Facts{{1, 1, 3, 2}, {2, 1, 3, 2}}), 9);  // all a_i == 2.

    // the same three, but asserting the oracle agrees with the literals too.
    CPF_EQ(reference(2, Facts{}), 2);
    CPF_EQ(reference(3, Facts{{1, 1, 3, 2}}), 5);
    CPF_EQ(reference(3, Facts{{1, 1, 3, 2}, {2, 1, 3, 2}}), 9);

    // --- edges the model has to nail. ---
    // q = 0: no facts, so spread n items over n values, one each -> cost n.
    CPF_EQ(solve(1, Facts{}), 1);
    CPF_EQ(solve(5, Facts{}), 5);
    CPF_EQ(solve(50, Facts{}), 50);

    // contradictory facts: a_i >= 2 and a_i <= 1 on the same range -> -1.
    CPF_EQ(solve(3, Facts{{1, 1, 3, 2}, {2, 1, 3, 1}}), -1);
    // partial contradiction: only position 2 is squeezed empty -> still -1.
    CPF_EQ(solve(4, Facts{{1, 2, 2, 3}, {2, 2, 2, 1}}), -1);

    // forced constant: every a_i pinned to one value v -> cnt(v) = n -> cost n^2.
    CPF_EQ(solve(4, Facts{{1, 1, 4, 3}, {2, 1, 4, 3}}), 16);
    CPF_EQ(solve(7, Facts{{1, 1, 7, 5}, {2, 1, 7, 5}}), 49);

    // n = 1: one item, one value. lower bound 2 on a size-1 array is impossible.
    CPF_EQ(solve(1, Facts{{1, 1, 1, 1}}), 1);
    CPF_EQ(solve(1, Facts{{1, 1, 1, 2}}), -1);  // a_1 >= 2 but v in [1,1].

    // half the array shoved into a band, the other half free -- the flow has to
    // balance the band without starving the free half.
    CPF_EQ(solve(4, Facts{{1, 1, 2, 2}, {2, 1, 2, 3}}),
           reference(4, Facts{{1, 1, 2, 2}, {2, 1, 2, 3}}));

    // --- thousands of tiny random cases vs the brute oracle. random facts drift
    //     into empty intervals on their own, so -1 gets exercised for free. ---
    auto sol = [](const Input& in) { return solve(in.n, in.facts); };
    auto ref = [](const Input& in) { return reference(in.n, in.facts); };

    // batch A: small n, up to a handful of facts. broadest coverage.
    {
        auto gen = [](cpf::Rng& r) {
            int n = static_cast<int>(r.in_range(1, 5));
            int q = static_cast<int>(r.in_range(0, 6));
            Facts f(static_cast<std::size_t>(q));
            for (auto& x : f) {
                x.t = static_cast<int>(r.in_range(1, 2));
                int a = static_cast<int>(r.in_range(1, n));
                int b = static_cast<int>(r.in_range(1, n));
                x.l = std::min(a, b);
                x.r = std::max(a, b);
                x.v = static_cast<int>(r.in_range(1, n));
            }
            return Input{n, f};
        };
        bool ok = cpf::differential(6000, 0xA1B05EEDULL, gen, sol, ref);
        CPF_CHECK(ok);
        std::printf("batch A (n<=5): 6000 cases\n");
    }

    // batch B: heavier fact load so infeasible -1 and tight forced values
    // dominate. tighter ranges keep the brute cheap.
    {
        auto gen = [](cpf::Rng& r) {
            int n = static_cast<int>(r.in_range(1, 6));
            int q = static_cast<int>(r.in_range(1, 12));
            Facts f(static_cast<std::size_t>(q));
            for (auto& x : f) {
                x.t = static_cast<int>(r.in_range(1, 2));
                int a = static_cast<int>(r.in_range(1, n));
                int b = static_cast<int>(r.in_range(1, n));
                x.l = std::min(a, b);
                x.r = std::max(a, b);
                x.v = static_cast<int>(r.in_range(1, n));
            }
            return Input{n, f};
        };
        bool ok = cpf::differential(6000, 0xDEADBEEFULL, gen, sol, ref);
        CPF_CHECK(ok);
        std::printf("batch B (n<=6): 6000 cases\n");
    }

    // batch C: push n to 7 (the oracle's exponential cap) with few facts, so the
    // widest intervals still get walked. fewer cases -- 7^7 leaves is real work.
    {
        auto gen = [](cpf::Rng& r) {
            int n = static_cast<int>(r.in_range(6, 7));
            int q = static_cast<int>(r.in_range(0, 4));
            Facts f(static_cast<std::size_t>(q));
            for (auto& x : f) {
                x.t = static_cast<int>(r.in_range(1, 2));
                int a = static_cast<int>(r.in_range(1, n));
                int b = static_cast<int>(r.in_range(1, n));
                x.l = std::min(a, b);
                x.r = std::max(a, b);
                x.v = static_cast<int>(r.in_range(1, n));
            }
            return Input{n, f};
        };
        bool ok = cpf::differential(400, 0xC0FFEEULL, gen, sol, ref);
        CPF_CHECK(ok);
        std::printf("batch C (n<=7): 400 cases\n");
    }

    return cpf::report();
}
