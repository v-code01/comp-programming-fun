#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using ll = std::int64_t;
using Edge = std::tuple<int, int, int, ll>;

struct Instance {
    int n = 0;
    std::vector<Edge> edges;
};

// a random weighted colored graph. FEW colors on purpose -- the partition cap
// has to bite. MIXED-SIGN weights on purpose -- so the empty set and partial
// sets win often, and the max-profit path selection earns its keep. u,v drawn
// independently, so self-loops and parallel edges appear on their own.
Instance random_graph(cpf::Rng& rng, int nlo, int nhi, int mlo, int mhi,
                      int ncolors, ll wlo, ll whi) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(nlo, nhi));
    const int m = static_cast<int>(rng.in_range(mlo, mhi));
    for (int e = 0; e < m; ++e) {
        const int u = static_cast<int>(rng.in_range(1, in.n));
        const int v = static_cast<int>(rng.in_range(1, in.n));
        const int c = static_cast<int>(rng.in_range(1, ncolors));
        const ll w = rng.in_range(wlo, whi);
        in.edges.emplace_back(u, v, c, w);
    }
    return in;
}

ll sol(const Instance& in) {
    return p0053::max_weight_colorful_forest(in.n, in.edges);
}
ll ref(const Instance& in) {
    return p0053::max_weight_colorful_forest_ref(in.n, in.edges);
}

}  // namespace

int main() {
    using namespace p0053;

    // ---- hand shapes. every expected number is the true answer AND matches the
    // 2^m oracle, so nothing here rests on my arithmetic alone. ----

    // a path 1-2-3-4, three distinct colors, all positive. already a forest,
    // already colorful -> take all of it. 5+3+7 = 15.
    {
        std::vector<Edge> e = {{1, 2, 1, 5}, {2, 3, 2, 3}, {3, 4, 3, 7}};
        CPF_EQ(max_weight_colorful_forest(4, e), 15);
        CPF_EQ(max_weight_colorful_forest_ref(4, e), 15);
    }

    // every weight negative -> touching any edge only loses weight. empty wins. 0.
    {
        std::vector<Edge> e = {{1, 2, 1, -5}, {2, 3, 2, -3}, {3, 4, 3, -7}};
        CPF_EQ(max_weight_colorful_forest(4, e), 0);
        CPF_EQ(max_weight_colorful_forest_ref(4, e), 0);
    }

    // the color cap forces dropping a heavy edge. edges 0 and 1 share color 1;
    // you keep only the heavier (10, not 4). edge 2 is color 2. best colorful
    // forest = {1-2 c1 w10, 3-4 c2 w8} = 18, disjoint so it IS a forest.
    {
        std::vector<Edge> e = {{1, 2, 1, 10}, {2, 3, 1, 4}, {3, 4, 2, 8}};
        CPF_EQ(max_weight_colorful_forest(4, e), 18);
        CPF_EQ(max_weight_colorful_forest_ref(4, e), 18);
    }

    // a SMALLER set beats the max-cardinality one. the lone +100 edge is a size-1
    // set worth 100. any size-2 colorful forest here drags in a -1, topping out
    // at 99. so the peak of W(k) sits at k=1, not the biggest k. -> 100.
    {
        std::vector<Edge> e = {{1, 2, 1, 100}, {2, 3, 2, -1}, {1, 3, 3, -1}};
        CPF_EQ(max_weight_colorful_forest(3, e), 100);
        CPF_EQ(max_weight_colorful_forest_ref(3, e), 100);
    }

    // both matroids bite at once, mixed signs. the oracle settles the number.
    {
        std::vector<Edge> e = {{1, 2, 1, 6}, {2, 3, 2, 5}, {3, 1, 1, 9},
                               {3, 4, 2, -4}, {4, 1, 3, 3}};
        CPF_EQ(max_weight_colorful_forest(4, e), max_weight_colorful_forest_ref(4, e));
        // color1: pick the heavier of {1-2 w6, 3-1 w9} -> 3-1 w9. color2: heavier
        // of {2-3 w5, 3-4 w-4} -> 2-3 w5. color3: 4-1 w3. {3-1,2-3,4-1} is a tree
        // on {1,2,3,4} -> 9+5+3 = 17.
        CPF_CHECK(max_weight_colorful_forest(4, e) == 17);
    }

    // ---- edges. ----

    // no edges at all -> nothing to pick, 0.
    {
        std::vector<Edge> e = {};
        CPF_EQ(max_weight_colorful_forest(5, e), 0);
        CPF_EQ(max_weight_colorful_forest_ref(5, e), 0);
    }

    // every edge the same color -> at most one survives, so the answer is the
    // single heaviest weight, here 9. negatives among them never get picked.
    {
        std::vector<Edge> e = {{1, 2, 7, 3}, {2, 3, 7, 9}, {3, 4, 7, -2}};
        CPF_EQ(max_weight_colorful_forest(4, e), 9);
        CPF_EQ(max_weight_colorful_forest_ref(4, e), 9);
    }

    // all one color AND all negative -> the single heaviest is still a loss, so
    // the empty set wins. 0.
    {
        std::vector<Edge> e = {{1, 2, 3, -3}, {2, 3, 3, -9}, {3, 4, 3, -2}};
        CPF_EQ(max_weight_colorful_forest(4, e), 0);
        CPF_EQ(max_weight_colorful_forest_ref(4, e), 0);
    }

    // a triangle, distinct colors, all positive. a forest takes any 2 of the 3
    // edges; take the two heaviest. 5+7 = 12 (drop the 3-weight edge).
    {
        std::vector<Edge> e = {{1, 2, 1, 5}, {2, 3, 2, 3}, {1, 3, 3, 7}};
        CPF_EQ(max_weight_colorful_forest(3, e), 12);
        CPF_EQ(max_weight_colorful_forest_ref(3, e), 12);
    }

    // a lone self-loop is never part of a forest -> 0, however heavy.
    {
        std::vector<Edge> e = {{2, 2, 1, 1000000}};
        CPF_EQ(max_weight_colorful_forest(3, e), 0);
        CPF_EQ(max_weight_colorful_forest_ref(3, e), 0);
    }

    // parallel edges of DIFFERENT colors between the same pair. only one can be
    // in a forest (the second closes a 2-cycle) -> take the heaviest, 8.
    {
        std::vector<Edge> e = {{1, 2, 1, 5}, {1, 2, 2, 8}, {1, 2, 3, 2}};
        CPF_EQ(max_weight_colorful_forest(2, e), 8);
        CPF_EQ(max_weight_colorful_forest_ref(2, e), 8);
    }

    // int64 reach: a star on 101 vertices, 100 distinct colors, every edge +1e6.
    // it's already a colorful forest, so the answer is the full sum 1e8 -- past
    // where a 32-bit accumulator would still be fine, but the point is the width
    // is exercised end to end. m=100 is far beyond the 2^m oracle, so this one is
    // checked against the closed-form sum, not by exhaustion.
    {
        std::vector<Edge> e;
        for (int k = 0; k < 100; ++k)
            e.emplace_back(1, k + 2, k + 1, static_cast<ll>(1000000));
        CPF_EQ(max_weight_colorful_forest(101, e), static_cast<ll>(100000000));
    }

    // ---- randomized differential vs the 2^m oracle. thousands of small weighted
    // colored graphs. small m keeps 2^m cheap; few colors and mixed-sign weights
    // keep both matroids and the empty/partial-set cases live. ----

    // sweep A: broad, 3 colors, small mixed weights -> ties, cap binds, cycles.
    {
        const int kCases = 30000;
        bool ok = cpf::differential(
            kCases, 0x53A0u,
            [](cpf::Rng& r) { return random_graph(r, 1, 7, 0, 12, 3, -10, 10); },
            sol, ref);
        CPF_CHECK(ok);
        std::printf("sweep A (n<=7, m<=12, 3 colors, w in [-10,10]): %d cases\n",
                    kCases);
    }

    // sweep B: only 2 colors -> the cap binds HARD, the answer rarely reaches a
    // spanning forest. denser edges so the graphic matroid rejects a lot. count
    // how often a smaller-than-max set (or the empty set) is the winner.
    {
        const int kCases = 20000;
        int diffs = 0, empty_wins = 0;
        cpf::Rng rng(0x53B0u);
        for (int t = 0; t < kCases; ++t) {
            Instance in = random_graph(rng, 2, 6, 0, 14, 2, -15, 15);
            const ll got = sol(in), want = ref(in);
            if (got != want) ++diffs;
            if (want == 0 && !in.edges.empty()) ++empty_wins;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep B (n<=6, m<=14, 2 colors, w in [-15,15]): %d cases, %d diffs, "
            "%d empty-set winners\n",
            kCases, diffs, empty_wins);
    }

    // sweep C: DENSE and small -> near-complete graphs, fundamental cycles
    // everywhere, hammering the exchange-graph tree-path scans. wider weights.
    {
        const int kCases = 15000;
        int diffs = 0, positive = 0;
        cpf::Rng rng(0x53C0u);
        for (int t = 0; t < kCases; ++t) {
            Instance in = random_graph(rng, 3, 5, 4, 12, 4, -8, 12);
            const ll got = sol(in), want = ref(in);
            if (got != want) ++diffs;
            if (want > 0) ++positive;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep C (n<=5 dense, m<=12, 4 colors, w in [-8,12]): %d cases, %d "
            "diffs, %d positive\n",
            kCases, diffs, positive);
    }

    // sweep D: single color -> the answer is the heaviest single edge or 0. a
    // sharp corner where the partition matroid alone decides. still diffed.
    {
        const int kCases = 10000;
        bool ok = cpf::differential(
            kCases, 0x53D0u,
            [](cpf::Rng& r) { return random_graph(r, 1, 8, 0, 13, 1, -12, 12); },
            sol, ref);
        CPF_CHECK(ok);
        std::printf("sweep D (n<=8, m<=13, 1 color, w in [-12,12]): %d cases\n",
                    kCases);
    }

    // sweep E: ALL-POSITIVE weights, more colors -> the cap loosens, larger sets
    // win, the augmenting loop runs deep. exercises the concave climb of W(k).
    {
        const int kCases = 15000;
        bool ok = cpf::differential(
            kCases, 0x53E0u,
            [](cpf::Rng& r) { return random_graph(r, 2, 7, 0, 13, 6, 1, 20); },
            sol, ref);
        CPF_CHECK(ok);
        std::printf("sweep E (n<=7, m<=13, 6 colors, w in [1,20]): %d cases\n",
                    kCases);
    }

    // sweep F: WIDE weights near the |w| bound -> the int64 accumulation is real.
    // small m so 2^m stays cheap while individual weights are big and signed.
    {
        const int kCases = 8000;
        int diffs = 0;
        cpf::Rng rng(0x53F0u);
        for (int t = 0; t < kCases; ++t) {
            Instance in = random_graph(rng, 1, 6, 0, 11, 4, -1000000, 1000000);
            const ll got = sol(in), want = ref(in);
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep F (n<=6, m<=11, 4 colors, w in [-1e6,1e6]): %d cases, %d "
            "diffs\n",
            kCases, diffs);
    }

    return cpf::report();
}
