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
using Edge = std::tuple<int, int, ll>;

struct Instance {
    int n = 0;
    std::vector<Edge> edges;
};

// a random general graph: pick n, then admit each unordered pair with a coin
// weighted by `density` (out of 100) and hang a random weight on it. cycles --
// including odd ones -- show up on their own the moment density climbs.
Instance random_graph(cpf::Rng& rng, int nlo, int nhi, ll wlo, ll whi,
                      int density) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(nlo, nhi));
    for (int i = 0; i < in.n; ++i)
        for (int j = i + 1; j < in.n; ++j)
            if (rng.in_range(1, 100) <= density)
                in.edges.emplace_back(i, j, rng.in_range(wlo, whi));
    return in;
}

// an ODD ring plus a scatter of chords -- built to make blossoms form and
// expand. the ring alone is a single odd cycle; the chords stack more of them.
Instance odd_cycle_graph(cpf::Rng& rng, ll wlo, ll whi) {
    Instance in;
    in.n = 2 * static_cast<int>(rng.in_range(1, 5)) + 1;  // 3, 5, 7, 9, or 11.
    for (int i = 0; i < in.n; ++i)
        in.edges.emplace_back(i, (i + 1) % in.n, rng.in_range(wlo, whi));
    const int chords = static_cast<int>(rng.in_range(0, in.n));
    for (int c = 0; c < chords; ++c) {
        int a = static_cast<int>(rng.in_range(0, in.n - 1));
        int b = static_cast<int>(rng.in_range(0, in.n - 1));
        if (a != b) in.edges.emplace_back(a, b, rng.in_range(wlo, whi));
    }
    return in;
}

ll sol(const Instance& in) { return p0032::max_weight_matching(in.n, in.edges); }
ll ref(const Instance& in) {
    return p0032::max_weight_matching_ref(in.n, in.edges);
}

}  // namespace

int main() {
    using namespace p0032;

    // ---- hand shapes. every "expected" is the true answer AND cross-checked
    // against the 2^n oracle, so no number here is just asserted from my head. ----

    // a triangle. you can take at most one edge -- any two share a vertex -- so
    // the heaviest single edge wins. 1, 2, 3 -> 3.
    {
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 2}, {0, 2, 3}};
        CPF_EQ(max_weight_matching(3, e), 3);
        CPF_EQ(max_weight_matching_ref(3, e), 3);
    }

    // a 4-cycle 0-1-2-3-0. the two matchings are the opposite-edge pairs; the
    // heavy pair {0-1, 2-3} = 5+5 beats {1-2, 3-0} = 1+1. -> 10.
    {
        std::vector<Edge> e = {{0, 1, 5}, {1, 2, 1}, {2, 3, 5}, {3, 0, 1}};
        CPF_EQ(max_weight_matching(4, e), 10);
        CPF_EQ(max_weight_matching_ref(4, e), 10);
    }

    // K4, all six edges. a perfect matching takes two disjoint edges; the best
    // disjoint pair here is {0-1 (10), 2-3 (9)} = 19, beating 8+7 and 6+5.
    {
        std::vector<Edge> e = {{0, 1, 10}, {2, 3, 9}, {0, 2, 8}, {1, 3, 7},
                               {0, 3, 6},  {1, 2, 5}};
        CPF_EQ(max_weight_matching(4, e), 19);
        CPF_EQ(max_weight_matching_ref(4, e), 19);
    }

    // a 5-cycle, unit weights -- the textbook blossom. an odd cycle of 5 admits
    // at most 2 disjoint edges. -> 2.
    {
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1},
                               {3, 4, 1}, {4, 0, 1}};
        CPF_EQ(max_weight_matching(5, e), 2);
        CPF_EQ(max_weight_matching_ref(5, e), 2);
    }

    // a 5-cycle with a heavy chord, so the blossom forms AND has to be resolved
    // in favor of the chord. ring weight 2, chord 0-2 weight 10: best is
    // {0-2 (10), 3-4 (2)} = 12.
    {
        std::vector<Edge> e = {{0, 1, 2}, {1, 2, 2}, {2, 3, 2}, {3, 4, 2},
                               {4, 0, 2}, {0, 2, 10}};
        CPF_EQ(max_weight_matching(5, e), 12);
        CPF_EQ(max_weight_matching_ref(5, e), 12);
    }

    // single vertex, no edge -> nothing to match, 0.
    {
        std::vector<Edge> e = {};
        CPF_EQ(max_weight_matching(1, e), 0);
        CPF_EQ(max_weight_matching_ref(1, e), 0);
    }

    // no edges on several vertices -> the empty matching, 0.
    {
        std::vector<Edge> e = {};
        CPF_EQ(max_weight_matching(5, e), 0);
        CPF_EQ(max_weight_matching_ref(5, e), 0);
    }

    // one heavy edge among idle vertices -> take it. n=6, only 2-5 present.
    {
        std::vector<Edge> e = {{2, 5, 1000000000LL}};
        CPF_EQ(max_weight_matching(6, e), 1000000000LL);
        CPF_EQ(max_weight_matching_ref(6, e), 1000000000LL);
    }

    // duplicate edges must dedup by max, not sum. 0-1 given as 3 and 7 -> 7.
    {
        std::vector<Edge> e = {{0, 1, 3}, {0, 1, 7}};
        CPF_EQ(max_weight_matching(2, e), 7);
        CPF_EQ(max_weight_matching_ref(2, e), 7);
    }

    // ---- randomized differential vs the 2^n subset-DP oracle. thousands of
    // small general graphs; small n keeps the oracle's 2^n cheap. ----

    // sweep A: broad, medium density -> mixed cycles, some blossoms.
    {
        const int kCases = 40000;
        bool ok = cpf::differential(
            kCases, 0x0032u,
            [](cpf::Rng& r) { return random_graph(r, 1, 10, 1, 40, 50); }, sol,
            ref);
        CPF_CHECK(ok);
        std::printf("sweep A (broad, 50%% density) vs oracle: %d cases, n<=10\n",
                    kCases);
    }

    // sweep B: DENSE -> near-complete graphs, blossoms form and expand often.
    {
        const int kCases = 25000;
        int nonzero = 0, diffs = 0;
        cpf::Rng rng(0xB105u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_graph(rng, 2, 9, 1, 25, 90);
            const ll got = sol(in), want = ref(in);
            if (want > 0) ++nonzero;
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep B (90%% dense) vs oracle: %d cases, n<=9, %d nonzero, %d "
            "diffs\n",
            kCases, nonzero, diffs);
    }

    // sweep C: odd rings + chords -> every case carries at least one odd cycle,
    // so the shrink/expand path is hit relentlessly.
    {
        const int kCases = 20000;
        int nonzero = 0, diffs = 0;
        cpf::Rng rng(0x0DDCu);
        for (int c = 0; c < kCases; ++c) {
            Instance in = odd_cycle_graph(rng, 1, 30);
            const ll got = sol(in), want = ref(in);
            if (want > 0) ++nonzero;
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "sweep C (odd rings + chords) vs oracle: %d cases, n<=11, %d "
            "nonzero, %d diffs\n",
            kCases, nonzero, diffs);
    }

    // sweep D: weights at the 1e9 scale so the int64 labels and duals get
    // stressed at real bounds, still diffed against the oracle on small n.
    {
        const int kCases = 5000;
        bool ok = cpf::differential(
            kCases, 0x1E9u,
            [](cpf::Rng& r) {
                return random_graph(r, 1, 8, 1, 1000000000LL, 70);
            },
            sol, ref);
        CPF_CHECK(ok);
        std::printf(
            "sweep D (weights up to 1e9) vs oracle: %d cases, n<=8\n", kCases);
    }

    return cpf::report();
}
