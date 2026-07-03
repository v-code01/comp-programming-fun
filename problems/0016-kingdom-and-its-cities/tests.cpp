#include <cstdint>
#include <cstdio>
#include <tuple>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// a random tree on n nodes: node i (2..n) picks a parent uniformly in 1..i-1.
// that's every labeled tree reachable, and it's the shape the real input can take.
std::vector<std::pair<int, int>> random_tree(int n, cpf::Rng& rng) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 2; i <= n; ++i) {
        const int p = static_cast<int>(rng.in_range(1, i - 1));
        edges.emplace_back(p, i);
    }
    return edges;
}

// a random non-empty important set: each node included with prob ~ half, but never
// empty. small k stresses the -1 path and single-junction cases; large k stresses
// the full compression.
std::vector<int> random_important(int n, cpf::Rng& rng) {
    std::vector<int> imp;
    for (int v = 1; v <= n; ++v) {
        if (rng.in_range(0, 1)) imp.push_back(v);
    }
    if (imp.empty()) imp.push_back(static_cast<int>(rng.in_range(1, n)));
    return imp;
}

// one query through a freshly built solver -- the differential's "solution" side.
int solver_answer(int n, const std::vector<std::pair<int, int>>& edges,
                  std::vector<int> imp) {
    p0016::Kingdom k;
    k.build(n, edges);
    return k.query(imp);
}

}  // namespace

int main() {
    using namespace p0016;

    // ---- the hand-checked shapes from the brief. verify each against the brute
    // oracle so the "expected" number is never just asserted. ----
    {
        // path 1-2-3. important {1,3} -> capture 2 -> answer 1.
        const std::vector<std::pair<int, int>> path3 = {{1, 2}, {2, 3}};
        CPF_EQ(solver_answer(3, path3, {1, 3}), 1);
        CPF_EQ(reference(3, path3, {1, 3}), 1);
        // important {1,2} are adjacent -> -1.
        CPF_EQ(solver_answer(3, path3, {1, 2}), -1);
        CPF_EQ(reference(3, path3, {1, 2}), -1);
    }

    // ---- a star: center 1, leaves 2,3,4. all leaves important -> cut the center,
    // answer 1. ----
    {
        const std::vector<std::pair<int, int>> star = {{1, 2}, {1, 3}, {1, 4}};
        CPF_EQ(solver_answer(4, star, {2, 3, 4}), 1);
        CPF_EQ(reference(4, star, {2, 3, 4}), 1);
        // two leaves important -> still just cut the center, answer 1.
        CPF_EQ(solver_answer(4, star, {2, 3}), 1);
        CPF_EQ(reference(4, star, {2, 3}), 1);
        // center + a leaf are adjacent -> -1.
        CPF_EQ(solver_answer(4, star, {1, 2}), -1);
    }

    // ---- single important city -> nothing to separate, answer 0. ----
    {
        const std::vector<std::pair<int, int>> path3 = {{1, 2}, {2, 3}};
        CPF_EQ(solver_answer(3, path3, {2}), 0);
        CPF_EQ(reference(3, path3, {2}), 0);
    }

    // ---- longer path 1-2-3-4-5, important {1,3,5} -> capture 2 and 4, answer 2.
    // two junctions, each forced, no single node severs both. ----
    {
        const std::vector<std::pair<int, int>> path5 = {
            {1, 2}, {2, 3}, {3, 4}, {4, 5}};
        CPF_EQ(solver_answer(5, path5, {1, 3, 5}), 2);
        CPF_EQ(reference(5, path5, {1, 3, 5}), 2);
    }

    // ---- deferral pays off: a broom. spine 1-2-3, then 3 forks to 4 and 5.
    // important {1,4,5}: paths 1..4 and 1..5 both run through 3. cut node 3 once
    // and both are severed -- answer 1, not 2. this is the >=2-junction save. ----
    {
        const std::vector<std::pair<int, int>> broom = {
            {1, 2}, {2, 3}, {3, 4}, {3, 5}};
        CPF_EQ(solver_answer(5, broom, {1, 4, 5}), 1);
        CPF_EQ(reference(5, broom, {1, 4, 5}), 1);
    }

    // ---- all leaves of a deeper tree important, mixed with an impossible case.
    // caterpillar: spine 1-2-3-4, leaves 5..8 hanging off 2 and 3. ----
    {
        const std::vector<std::pair<int, int>> cat = {
            {1, 2}, {2, 3}, {3, 4}, {2, 5}, {2, 6}, {3, 7}, {3, 8}};
        // leaves 5,6,7,8 important -> cut 2 (severs 5,6 from the rest) and 3
        // (severs 7,8): the oracle knows the true min.
        CPF_EQ(solver_answer(8, cat, {5, 6, 7, 8}), reference(8, cat, {5, 6, 7, 8}));
        // include the spine so 2 and one of its leaf-neighbors are both important
        // -> adjacent -> -1.
        CPF_EQ(solver_answer(8, cat, {2, 5}), -1);
        CPF_EQ(reference(8, cat, {2, 5}), -1);
    }

    // ---- one built solver answering a batch of queries: proves per-query scratch
    // is reset cleanly and never bleeds between queries. diff every answer. ----
    {
        cpf::Rng rng(0xB00Cu);
        const int n = 11;
        const auto edges = random_tree(n, rng);
        Kingdom k;
        k.build(n, edges);
        int diffs = 0;
        for (int t = 0; t < 4000; ++t) {
            auto imp = random_important(n, rng);
            auto imp_copy = imp;
            const int got = k.query(imp);              // reused solver.
            const int want = reference(n, edges, imp_copy);
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("reused-solver batch: 4000 queries on one tree, %d diffs\n",
                    diffs);
    }

    // ---- randomized differential vs the subset-brute oracle. small n so the
    // 2^(n-k) oracle stays cheap; thousands of trees x important sets, including
    // the -1 cases and multi-junction saves the random sets produce constantly. ----
    {
        const int kCases = 20000;
        bool ok = cpf::differential(
            kCases, 0x6613u,
            [](cpf::Rng& rng) {
                const int n = static_cast<int>(rng.in_range(1, 11));
                auto edges = random_tree(n, rng);
                auto imp = random_important(n, rng);
                return std::make_tuple(n, edges, imp);
            },
            [](const auto& in) {
                auto imp = std::get<2>(in);  // solver sorts in place -- give it a copy.
                return solver_answer(std::get<0>(in), std::get<1>(in), imp);
            },
            [](const auto& in) {
                return reference(std::get<0>(in), std::get<1>(in), std::get<2>(in));
            });
        CPF_CHECK(ok);
        std::printf("randomized differential vs subset-brute: %d cases, n<=11\n",
                    kCases);
    }

    // ---- a second differential skewed tiny (n<=7) so -1 and dense important sets
    // dominate -- hammers the adjacency detection specifically. ----
    {
        const int kCases = 20000;
        int minus_ones = 0;
        int diffs = 0;
        cpf::Rng rng(0x1CE5u);
        for (int t = 0; t < kCases; ++t) {
            const int n = static_cast<int>(rng.in_range(2, 7));
            const auto edges = random_tree(n, rng);
            auto imp = random_important(n, rng);
            auto imp_copy = imp;
            const int got = solver_answer(n, edges, imp);
            const int want = reference(n, edges, imp_copy);
            if (want == -1) ++minus_ones;
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("tiny differential n<=7: %d cases, %d were -1, %d diffs\n",
                    kCases, minus_ones, diffs);
    }

    return cpf::report();
}
