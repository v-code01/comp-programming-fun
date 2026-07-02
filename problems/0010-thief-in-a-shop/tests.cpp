#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0010::solve;
using Vec = std::vector<int>;

// what the differential feeds both solvers: k and the price list.
struct Input {
    int k;
    Vec a;
};

int main() {
    // --- the three statement examples, whole lists compared. ---
    CPF_EQ(solve(2, Vec{1, 2, 3}), (Vec{2, 3, 4, 5, 6}));
    CPF_EQ(solve(5, Vec{1, 1, 1, 1, 1}), (Vec{5}));
    CPF_EQ(solve(3, Vec{3, 5, 11}),
           (Vec{9, 11, 13, 15, 17, 19, 21, 25, 27, 33}));

    // --- edges the shift-and-pad argument has to nail. ---
    // k = 1: exactly the distinct prices, sorted -- one pick is one price.
    CPF_EQ(solve(1, Vec{5, 2, 8, 2}), (Vec{2, 5, 8}));
    // every price equal: only k*a is reachable, the pad is the whole multiset.
    CPF_EQ(solve(4, Vec{7, 7, 7}), (Vec{28}));
    // n = 1: one kind, one total.
    CPF_EQ(solve(6, Vec{3}), (Vec{18}));
    // several kinds, single distinct value -- same as all-equal.
    CPF_EQ(solve(3, Vec{2, 2, 2, 2}), (Vec{6}));
    // min price is 1 (b has a real zero already) and a gap the coins can't fill:
    // coin is 4 (=5-1), so with k=2 the extras are {0,4,8} -> totals {2,6,10}.
    CPF_EQ(solve(2, Vec{1, 5}), (Vec{2, 6, 10}));

    // --- thousands of small random cases vs the exactly-k set oracle. ---
    // tiny n, k, costs -- the oracle's array stays trivial and the full sorted
    // lists must match. any disagreement prints its seed and stops.
    {
        auto gen = [](cpf::Rng& r) {
            int n = static_cast<int>(r.in_range(1, 6));
            int k = static_cast<int>(r.in_range(1, 8));
            Vec a(static_cast<std::size_t>(n));
            for (auto& x : a) x = static_cast<int>(r.in_range(1, 10));
            return Input{k, a};
        };
        auto sol = [](const Input& in) { return solve(in.k, in.a); };
        auto ref = [](const Input& in) { return p0010::reference(in.k, in.a); };
        bool ok = cpf::differential(6000, 0xB16B00B5ULL, gen, sol, ref);
        CPF_CHECK(ok);
        std::printf("small differential: 6000 cases\n");
    }

    // --- a second batch, skewed the other way: more kinds, bigger k, tight cost
    // range so collisions and gaps both show up often. ---
    {
        auto gen = [](cpf::Rng& r) {
            int n = static_cast<int>(r.in_range(1, 8));
            int k = static_cast<int>(r.in_range(1, 10));
            Vec a(static_cast<std::size_t>(n));
            for (auto& x : a) x = static_cast<int>(r.in_range(1, 5));
            return Input{k, a};
        };
        auto sol = [](const Input& in) { return solve(in.k, in.a); };
        auto ref = [](const Input& in) { return p0010::reference(in.k, in.a); };
        bool ok = cpf::differential(6000, 0xDEADC0DEULL, gen, sol, ref);
        CPF_CHECK(ok);
        std::printf("dense differential: 6000 cases\n");
    }

    return cpf::report();
}
