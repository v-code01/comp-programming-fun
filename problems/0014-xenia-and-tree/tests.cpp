#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0014::Bfs;
using p0014::NearestRed;

// one operation in the statement's alphabet, vertices already 0-based.
//   type 1: paint v red
//   type 2: query nearest-red distance from v
struct Op {
    int type;
    int v;
};

// replay a program on any structure with the build/paint/query API, collecting
// the type-2 answers in order. templated so the centroid tree and the oracle run
// the identical path -- build() paints vertex 0 for both.
template <typename DS>
std::vector<int> run(int n, const std::vector<std::pair<int, int>>& edges,
                     const std::vector<Op>& ops) {
    DS ds;
    ds.build(n, edges);
    std::vector<int> out;
    for (const Op& o : ops) {
        if (o.type == 1)
            ds.paint(o.v);
        else
            out.push_back(ds.query(o.v));
    }
    return out;
}

// centroid tree vs oracle on one fixed program -- both must emit the same
// type-2 answer stream.
void expect_agree(int n, const std::vector<std::pair<int, int>>& edges,
                  const std::vector<Op>& ops) {
    CPF_CHECK(run<NearestRed>(n, edges, ops) == run<Bfs>(n, edges, ops));
}

// a path 0-1-2-...-(n-1).
std::vector<std::pair<int, int>> path_edges(int n) {
    std::vector<std::pair<int, int>> e;
    for (int i = 1; i < n; ++i) e.emplace_back(i - 1, i);
    return e;
}

// a star: vertex 0 at the center, every other vertex a leaf off it.
std::vector<std::pair<int, int>> star_edges(int n) {
    std::vector<std::pair<int, int>> e;
    for (int i = 1; i < n; ++i) e.emplace_back(0, i);
    return e;
}

// a uniformly random tree: each vertex i in 1..n-1 picks a parent in [0, i-1].
// connected by construction, 0-based, no LCA/relabel games.
std::vector<std::pair<int, int>> random_tree(int n, cpf::Rng& rng) {
    std::vector<std::pair<int, int>> e;
    for (int i = 1; i < n; ++i) {
        e.emplace_back(static_cast<int>(rng.in_range(0, i - 1)), i);
    }
    return e;
}

}  // namespace

int main() {
    // ---- the statement example ------------------------------------------------
    // n=5, edges (1,2)(2,3)(2,4)(4,5); queries 2 1 -> 0, 2 5 -> 3, 1 2, 2 5 -> 2.
    // 1-based -> 0-based: edges (0,1)(1,2)(1,3)(3,4); vertex 1 -> 0, vertex 5 -> 4.
    {
        const int n = 5;
        std::vector<std::pair<int, int>> edges{{0, 1}, {1, 2}, {1, 3}, {3, 4}};
        std::vector<Op> ops{
            {2, 0},  // 2 1 -> 0 (vertex 1 is red from the start)
            {2, 4},  // 2 5 -> 3
            {1, 1},  // 1 2 (paint vertex 2)
            {2, 4},  // 2 5 -> 2
        };
        std::vector<int> got = run<NearestRed>(n, edges, ops);
        std::vector<int> want{0, 3, 2};
        CPF_CHECK(got == want);
        expect_agree(n, edges, ops);
    }

    // ---- enumerated edges -----------------------------------------------------

    // n=1: the lone vertex is vertex 1, already red. distance to itself is 0.
    {
        const int n = 1;
        std::vector<std::pair<int, int>> edges;
        std::vector<Op> ops{{2, 0}, {1, 0}, {2, 0}};  // 0, then still 0.
        CPF_CHECK(run<NearestRed>(n, edges, ops) == (std::vector<int>{0, 0}));
        expect_agree(n, edges, ops);
    }

    // n=2, single edge. query the far end (3), paint it, re-query (0).
    {
        const int n = 2;
        std::vector<std::pair<int, int>> edges{{0, 1}};
        std::vector<Op> ops{{2, 1}, {1, 1}, {2, 1}, {2, 0}};
        CPF_CHECK(run<NearestRed>(n, edges, ops) == (std::vector<int>{1, 0, 0}));
        expect_agree(n, edges, ops);
    }

    // a path -- the worst case for centroid-tree depth. query the far tip before
    // and after painting the midpoint, so the answer must halve.
    {
        const int n = 9;  // 0..8, midpoint 4.
        auto edges = path_edges(n);
        std::vector<Op> ops{
            {2, 8},  // far tip, only vertex 0 red -> 8
            {1, 4},  // paint the midpoint
            {2, 8},  // nearest red is now 4 -> 4
            {1, 8},  // paint the tip itself
            {2, 8},  // -> 0
            {2, 6},  // between 4 and 8, both red -> 2
        };
        CPF_CHECK(run<NearestRed>(n, edges, ops) ==
                  (std::vector<int>{8, 4, 0, 2}));
        expect_agree(n, edges, ops);
    }

    // a star -- every leaf is distance 1 from the center (vertex 0, red). painting
    // one leaf leaves the others at 2 (leaf -> center -> leaf).
    {
        const int n = 6;
        auto edges = star_edges(n);
        std::vector<Op> ops{
            {2, 3},  // leaf -> center is red -> 1
            {1, 3},  // paint leaf 3
            {2, 3},  // itself red -> 0
            {2, 5},  // another leaf: nearest red is center -> 1
        };
        CPF_CHECK(run<NearestRed>(n, edges, ops) == (std::vector<int>{1, 0, 1}));
        expect_agree(n, edges, ops);
    }

    // repeated paints of the same vertex must be idempotent -- the second paint
    // relaxes against distances already recorded and changes nothing.
    {
        const int n = 7;
        auto edges = path_edges(n);
        std::vector<Op> ops{
            {1, 5}, {1, 5}, {1, 5},  // paint 5 three times
            {2, 6},                  // -> 1
            {2, 0},                  // vertex 0 still red -> 0
        };
        CPF_CHECK(run<NearestRed>(n, edges, ops) == (std::vector<int>{1, 0}));
        expect_agree(n, edges, ops);
    }

    // query with no paints at all: only vertex 0 is red, answers are plain
    // distances from 0. exercised against the oracle on a path.
    {
        const int n = 6;
        auto edges = path_edges(n);
        std::vector<Op> ops{{2, 0}, {2, 1}, {2, 3}, {2, 5}};
        CPF_CHECK(run<NearestRed>(n, edges, ops) ==
                  (std::vector<int>{0, 1, 3, 5}));
        expect_agree(n, edges, ops);
    }

    // ---- randomized differential vs the BFS oracle ----------------------------
    // thousands of independent programs on random trees. small n so the O(n)-per-
    // query oracle stays cheap; paint-heavy and query-heavy mixes both appear so
    // best[] relaxations and the query min each get hammered. any disagreement
    // prints the seed that replays it.
    {
        constexpr int kProgs = 8000;
        cpf::Rng rng(0x0014ULL);
        int diffs = 0;
        std::uint64_t first_bad = 0;

        for (int s = 0; s < kProgs; ++s) {
            const int n = static_cast<int>(rng.in_range(1, 40));
            const auto edges = random_tree(n, rng);

            const int m = static_cast<int>(rng.in_range(1, 60));
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int q = 0; q < m; ++q) {
                // skew the mix per program: some are paint-heavy, some query-heavy.
                const int paint_pct = static_cast<int>(rng.in_range(20, 80));
                const bool paint = rng.in_range(0, 99) < paint_pct;
                const int v = static_cast<int>(rng.in_range(0, n - 1));
                ops.push_back({paint ? 1 : 2, v});
            }

            if (run<NearestRed>(n, edges, ops) != run<Bfs>(n, edges, ops)) {
                if (diffs == 0) first_bad = static_cast<std::uint64_t>(s);
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first differing program index: %llu\n",
                        static_cast<unsigned long long>(first_bad));
        }
        std::printf("differential: %d programs, %d diffs\n", kProgs, diffs);
    }

    // ---- a few larger random trees, query-dense, still vs the oracle ----------
    // bigger n exercises deeper centroid trees -- more ancestors per vertex, more
    // best[] entries interacting. the oracle stays O(n) per query so keep counts
    // modest but push n up.
    {
        constexpr int kProgs = 300;
        cpf::Rng rng(0xBEEFULL);
        int diffs = 0;
        for (int s = 0; s < kProgs; ++s) {
            const int n = static_cast<int>(rng.in_range(50, 400));
            const auto edges = random_tree(n, rng);
            const int m = 200;
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int q = 0; q < m; ++q) {
                const bool paint = rng.in_range(0, 2) == 0;  // ~1/3 paints.
                const int v = static_cast<int>(rng.in_range(0, n - 1));
                ops.push_back({paint ? 1 : 2, v});
            }
            if (run<NearestRed>(n, edges, ops) != run<Bfs>(n, edges, ops)) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("larger-tree differential: %d programs, %d diffs\n", kProgs,
                    diffs);
    }

    return cpf::report();
}
