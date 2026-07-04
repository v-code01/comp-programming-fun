#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0026::Offer;

struct Graph {
    int n;
    std::vector<Offer> offers;
    bool arb;  // constructed ground truth, margin-clear -- float can't flip it.
};

// a graph with no arbitrage, and the margin is not close. hand every node a
// potential, then price every edge strictly under its potential gap:
//   rate = exp(pot[u] - pot[v]) shrunk by exp(-delta), delta in [0.05, 0.5].
// around any cycle the potentials telescope to 1, so the product is exp(-sum
// delta) <= exp(-0.05 * len) -- below 1 with room. rates stay in [1e-6, 1e6]
// because pot in [-2, 2] keeps the exponent inside [-4.5, 4].
Graph gen_safe(cpf::Rng& rng) {
    const int n = static_cast<int>(rng.in_range(2, 8));
    std::vector<double> pot(static_cast<std::size_t>(n) + 1, 0.0);
    for (int i = 1; i <= n; ++i)
        pot[static_cast<std::size_t>(i)] =
            static_cast<double>(rng.in_range(-200, 200)) / 100.0;

    Graph g{n, {}, false};
    const int m = static_cast<int>(rng.in_range(0, n * n));
    for (int e = 0; e < m; ++e) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        const double delta = static_cast<double>(rng.in_range(5, 50)) / 100.0;
        const double rate = std::exp(pot[static_cast<std::size_t>(u)] -
                                     pot[static_cast<std::size_t>(v)] - delta);
        g.offers.push_back({u, v, rate});
    }
    return g;
}

// a safe base with one loud cycle bolted on. the injected loop's product lands
// at 1.69 or 1.728 -- an arb no rounding noise can erase. adding a negative cycle
// to an arb-free graph makes the answer YES, full stop.
Graph gen_arb(cpf::Rng& rng) {
    Graph g = gen_safe(rng);
    g.arb = true;
    if (g.n >= 3 && rng.in_range(0, 1) == 0) {
        // 3-cycle 1->2->3->1, each rate 1.2 -> product 1.728.
        g.offers.push_back({1, 2, 1.2});
        g.offers.push_back({2, 3, 1.2});
        g.offers.push_back({3, 1, 1.2});
    } else {
        // 2-cycle 1->2->1, each rate 1.3 -> product 1.69.
        g.offers.push_back({1, 2, 1.3});
        g.offers.push_back({2, 1, 1.3});
    }
    return g;
}

}  // namespace

int main() {
    using namespace p0026;

    // ---- classic hand shapes. every verdict cross-checked against the
    // floyd-warshall oracle so no number is just asserted from my head. ----

    // three currencies, an obvious loop: 1->2 (2.0), 2->3 (1.0), 3->1 (0.8).
    // product 1.6 -- a dollar comes back $1.60. arbitrage.
    {
        std::vector<Offer> g = {{1, 2, 2.0}, {2, 3, 1.0}, {3, 1, 0.8}};
        CPF_EQ(has_arbitrage(3, g), true);
        CPF_EQ(has_arbitrage_reference(3, g), true);
    }
    // three currencies priced consistently and shrinking: every leg 0.9 both
    // ways. the tightest loop is 0.9*0.9 = 0.81 < 1. no arbitrage.
    {
        std::vector<Offer> g = {{1, 2, 0.9}, {2, 1, 0.9}, {2, 3, 0.9},
                                {3, 2, 0.9}, {1, 3, 0.9}, {3, 1, 0.9}};
        CPF_EQ(has_arbitrage(3, g), false);
        CPF_EQ(has_arbitrage_reference(3, g), false);
    }

    // ---- edges the algorithm has to nail. ----

    // no edges at all -> nothing to loop through. NO.
    {
        std::vector<Offer> g = {};
        CPF_EQ(has_arbitrage(1, g), false);
        CPF_EQ(has_arbitrage_reference(1, g), false);
        CPF_EQ(has_arbitrage(300, g), false);
        CPF_EQ(has_arbitrage_reference(300, g), false);
    }
    // a self-loop that mints money: 1->1 at rate 2.0. a length-1 negative cycle,
    // caught by the detection sweep alone. YES.
    {
        std::vector<Offer> g = {{1, 1, 2.0}};
        CPF_EQ(has_arbitrage(1, g), true);
        CPF_EQ(has_arbitrage_reference(1, g), true);
    }
    // a self-loop that bleeds money: 1->1 at rate 0.5. product 0.5 < 1. NO.
    {
        std::vector<Offer> g = {{1, 1, 0.5}};
        CPF_EQ(has_arbitrage(1, g), false);
        CPF_EQ(has_arbitrage_reference(1, g), false);
    }
    // two nodes, u->v->u with product > 1: 1->2 (1.5), 2->1 (1.0) -> 1.5. YES.
    {
        std::vector<Offer> g = {{1, 2, 1.5}, {2, 1, 1.0}};
        CPF_EQ(has_arbitrage(2, g), true);
        CPF_EQ(has_arbitrage_reference(2, g), true);
    }
    // two nodes, product under 1: 1->2 (1.0), 2->1 (0.5) -> 0.5. NO.
    {
        std::vector<Offer> g = {{1, 2, 1.0}, {2, 1, 0.5}};
        CPF_EQ(has_arbitrage(2, g), false);
        CPF_EQ(has_arbitrage_reference(2, g), false);
    }
    // a one-way chain with no loop back: 1->2->3, no return edge. no cycle at
    // all, so no arb however good the rates look. NO.
    {
        std::vector<Offer> g = {{1, 2, 5.0}, {2, 3, 5.0}};
        CPF_EQ(has_arbitrage(3, g), false);
        CPF_EQ(has_arbitrage_reference(3, g), false);
    }
    // rates at the bound: 1->2 at 1e6, 2->1 at 1e-6 -> product exactly balanced
    // under the shrink, no arb. pushes the -log magnitude to its extreme.
    {
        std::vector<Offer> g = {{1, 2, 1000000.0}, {2, 1, 0.000001}};
        CPF_EQ(has_arbitrage(2, g), false);
        CPF_EQ(has_arbitrage_reference(2, g), false);
    }

    // ---- randomized differential vs the floyd-warshall oracle. thousands of
    // margin-clear small graphs, half arb-free and half with a loud injected
    // loop, so the YES/NO ground truth can't be flipped by float error and both
    // algorithms land on it. ----
    {
        const int kCases = 40000;
        cpf::Rng rng(0x2626u);
        int diffs = 0;         // solution vs oracle.
        int mislabels = 0;     // either one vs the constructed truth.
        int yes = 0, no = 0;
        int first_diff_seed = -1;
        for (int c = 0; c < kCases; ++c) {
            const Graph g = (c & 1) ? gen_arb(rng) : gen_safe(rng);
            const bool sol = has_arbitrage(g.n, g.offers);
            const bool ref = has_arbitrage_reference(g.n, g.offers);
            if (sol != ref && first_diff_seed < 0) first_diff_seed = c;
            if (sol != ref) ++diffs;
            if (sol != g.arb || ref != g.arb) ++mislabels;
            g.arb ? ++yes : ++no;
        }
        CPF_CHECK(diffs == 0);
        CPF_CHECK(mislabels == 0);
        if (first_diff_seed >= 0)
            std::printf("  first diff at case %d\n", first_diff_seed);
        std::printf("randomized differential vs floyd-warshall: %d cases "
                    "(%d arb, %d safe), %d diffs, %d mislabels\n",
                    kCases, yes, no, diffs, mislabels);
    }

    return cpf::report();
}
