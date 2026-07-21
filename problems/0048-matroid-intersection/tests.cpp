#include <cstdio>
#include <tuple>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Edge = std::tuple<int, int, int>;

struct Instance {
    int n = 0;
    std::vector<Edge> edges;
};

// a random colored graph. FEW colors on purpose -- the partition cap has to
// actually bind, so the answer lands below the spanning-forest size and the
// exchange graph does real work. u,v drawn independently, so self-loops and
// parallel edges show up on their own, exactly what the constraints allow.
Instance random_graph(cpf::Rng& rng, int nlo, int nhi, int mlo, int mhi,
                      int ncolors) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(nlo, nhi));
    const int m = static_cast<int>(rng.in_range(mlo, mhi));
    for (int e = 0; e < m; ++e) {
        const int u = static_cast<int>(rng.in_range(1, in.n));
        const int v = static_cast<int>(rng.in_range(1, in.n));
        const int c = static_cast<int>(rng.in_range(1, ncolors));
        in.edges.emplace_back(u, v, c);
    }
    return in;
}

int sol(const Instance& in) { return p0048::max_colorful_forest(in.n, in.edges); }
int ref(const Instance& in) {
    return p0048::max_colorful_forest_ref(in.n, in.edges);
}

}  // namespace

int main() {
    using namespace p0048;

    // ---- hand shapes. every expected number is the true answer AND matches the
    // 2^m oracle, so nothing here is asserted from my head alone. ----

    // a triangle, all one color. the forest would take 2 of the 3 edges, but one
    // color caps it at a single edge. -> 1.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 3, 1}, {1, 3, 1}};
        CPF_EQ(max_colorful_forest(3, e), 1);
        CPF_EQ(max_colorful_forest_ref(3, e), 1);
    }

    // a path 1-2-3-4 with three distinct colors. it's already a forest and
    // already colorful -- take all of it. n-1 = 3.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 3, 2}, {3, 4, 3}};
        CPF_EQ(max_colorful_forest(4, e), 3);
        CPF_EQ(max_colorful_forest_ref(4, e), 3);
    }

    // the color cap binds BELOW the spanning forest. n=4 admits a 3-edge forest,
    // but only two colors exist across the graph, so at most 2 edges survive --
    // {1-2 color 1, 3-4 color 2} is a disjoint colorful forest. -> 2, not 3.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 3, 1}, {3, 4, 2}, {1, 4, 2}};
        CPF_EQ(max_colorful_forest(4, e), 2);
        CPF_EQ(max_colorful_forest_ref(4, e), 2);
    }

    // ---- edges. ----

    // no edges at all -> nothing to pick, 0.
    {
        std::vector<Edge> e = {};
        CPF_EQ(max_colorful_forest(5, e), 0);
        CPF_EQ(max_colorful_forest_ref(5, e), 0);
    }

    // every edge the same color -> at most one edge, whatever the shape. 1.
    {
        std::vector<Edge> e = {{1, 2, 7}, {2, 3, 7}, {3, 4, 7}, {4, 5, 7}};
        CPF_EQ(max_colorful_forest(5, e), 1);
        CPF_EQ(max_colorful_forest_ref(5, e), 1);
    }

    // a forest that is already colorful returns its own size. a star on 5
    // vertices, four distinct colors -> 4.
    {
        std::vector<Edge> e = {{1, 2, 1}, {1, 3, 2}, {1, 4, 3}, {1, 5, 4}};
        CPF_EQ(max_colorful_forest(5, e), 4);
        CPF_EQ(max_colorful_forest_ref(5, e), 4);
    }

    // a lone self-loop is never part of a forest -> 0.
    {
        std::vector<Edge> e = {{2, 2, 1}};
        CPF_EQ(max_colorful_forest(3, e), 0);
        CPF_EQ(max_colorful_forest_ref(3, e), 0);
    }

    // parallel edges of DIFFERENT colors between the same pair -- still only one
    // can be in a forest (the second closes a 2-cycle). -> 1.
    {
        std::vector<Edge> e = {{1, 2, 1}, {1, 2, 2}, {1, 2, 3}};
        CPF_EQ(max_colorful_forest(2, e), 1);
        CPF_EQ(max_colorful_forest_ref(2, e), 1);
    }

    // a case where both constraints bite at once. two triangles sharing colors:
    // 4-cycle-ish with a heavy color reuse. the oracle settles the number.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 3, 2}, {3, 1, 1},
                               {3, 4, 2}, {4, 1, 3}};
        CPF_EQ(max_colorful_forest(4, e), max_colorful_forest_ref(4, e));
        CPF_CHECK(max_colorful_forest(4, e) == 3);  // {1-2 c1, 2-3 c2, 4-1 c3}.
    }

    // ---- randomized differential vs the 2^m oracle. thousands of small colored
    // graphs, few colors so the partition matroid actually binds and cycles force
    // graphic-matroid rejections. small m keeps the oracle's 2^m cheap. ----

    // sweep A: broad, 3 colors, moderate edge count -> cap binds, cycles common.
    {
        const int kCases = 30000;
        bool ok = cpf::differential(
            kCases, 0x0048u,
            [](cpf::Rng& r) { return random_graph(r, 1, 7, 0, 12, 3); }, sol, ref);
        CPF_CHECK(ok);
        std::printf("sweep A (n<=7, m<=12, 3 colors) vs oracle: %d cases\n",
                    kCases);
    }

    // sweep B: only 2 colors -> the partition cap binds HARD, answer rarely
    // reaches n-1. denser edges so the graphic matroid rejects a lot.
    {
        const int kCases = 15000;
        int diffs = 0, bound = 0;
        cpf::Rng rng(0xB2u);
        for (int t = 0; t < kCases; ++t) {
            Instance in = random_graph(rng, 2, 6, 0, 14, 2);
            const int got = sol(in), want = ref(in);
            if (got != want) ++diffs;
            if (want < in.n - 1) ++bound;  // cap or cycles held it below forest size.
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep B (n<=6, m<=14, 2 colors) vs oracle: %d cases, %d diffs, "
            "%d held below n-1\n",
            kCases, diffs, bound);
    }

    // sweep C: DENSE and small -> near-complete graphs where fundamental cycles
    // are everywhere, hammering the exchange-graph tree-path scans.
    {
        const int kCases = 12000;
        int diffs = 0, nonzero = 0;
        cpf::Rng rng(0xC3u);
        for (int t = 0; t < kCases; ++t) {
            Instance in = random_graph(rng, 3, 5, 4, 12, 4);
            const int got = sol(in), want = ref(in);
            if (got != want) ++diffs;
            if (want > 0) ++nonzero;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep C (n<=5 dense, m<=12, 4 colors) vs oracle: %d cases, %d diffs, "
            "%d nonzero\n",
            kCases, diffs, nonzero);
    }

    // sweep D: single color everywhere -> the answer is always 0 or 1, a sharp
    // corner where the partition matroid alone decides. still diffed.
    {
        const int kCases = 8000;
        bool ok = cpf::differential(
            kCases, 0xD4u,
            [](cpf::Rng& r) { return random_graph(r, 1, 8, 0, 13, 1); }, sol, ref);
        CPF_CHECK(ok);
        std::printf("sweep D (n<=8, m<=13, 1 color) vs oracle: %d cases\n",
                    kCases);
    }

    return cpf::report();
}
