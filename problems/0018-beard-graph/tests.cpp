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

using p0018::BeardHLD;
using p0018::BeardWalk;

// one op in the statement's alphabet, indices already 0-based.
//   type 1: paint edge e black
//   type 2: paint edge e white
//   type 3: query path a -> b   (e unused for the endpoints; a, b carry them)
struct Op {
    int type;
    int a;  // edge index (types 1,2) or first vertex (type 3).
    int b;  // second vertex (type 3 only).
};

// replay a program on any structure with the Beard API, collecting the type-3
// answers in order. templated so the HLD tree and the oracle run the identical
// path -- build() starts every edge black for both.
template <typename DS>
std::vector<int> run(int n, const std::vector<std::pair<int, int>>& edges,
                     const std::vector<Op>& ops) {
    DS ds;
    ds.build(n, edges);
    std::vector<int> out;
    for (const Op& o : ops) {
        if (o.type == 1)
            ds.paint_black(o.a);
        else if (o.type == 2)
            ds.paint_white(o.a);
        else
            out.push_back(ds.query(o.a, o.b));
    }
    return out;
}

// HLD vs oracle on one fixed program -- both must emit the same type-3 stream.
void expect_agree(int n, const std::vector<std::pair<int, int>>& edges,
                  const std::vector<Op>& ops) {
    CPF_CHECK(run<BeardHLD>(n, edges, ops) == run<BeardWalk>(n, edges, ops));
}

// a path 0-1-2-...-(n-1); edge i joins i and i+1.
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
// connected by construction, 0-based, edge index = i-1.
std::vector<std::pair<int, int>> random_tree(int n, cpf::Rng& rng) {
    std::vector<std::pair<int, int>> e;
    for (int i = 1; i < n; ++i) {
        e.emplace_back(static_cast<int>(rng.in_range(0, i - 1)), i);
    }
    return e;
}

}  // namespace

int main() {
    // ---- the statement example, brute-verified by hand ------------------------
    // n=3, path 1-2-3: edge1=(1,2), edge2=(2,3). queries, 1-based ->
    //   3 1 2 -> path {e1}, black -> 1
    //   3 1 3 -> path {e1,e2}, black -> 2
    //   3 2 3 -> path {e2}, black -> 1
    //   2 2   -> paint edge2 white
    //   3 1 2 -> {e1} black -> 1
    //   3 1 3 -> {e1,e2}, e2 white -> -1
    //   3 2 3 -> {e2} white -> -1
    {
        const int n = 3;
        std::vector<std::pair<int, int>> edges{{0, 1}, {1, 2}};
        std::vector<Op> ops{
            {3, 0, 1}, {3, 0, 2}, {3, 1, 2},  // -> 1, 2, 1
            {2, 1, 0},                          // paint edge index 1 (edge2) white
            {3, 0, 1}, {3, 0, 2}, {3, 1, 2},  // -> 1, -1, -1
        };
        std::vector<int> got = run<BeardHLD>(n, edges, ops);
        std::vector<int> want{1, 2, 1, 1, -1, -1};
        CPF_CHECK(got == want);
        expect_agree(n, edges, ops);
    }

    // ---- enumerated edges -----------------------------------------------------

    // n=1: a lone vertex, no edges. the only legal query is a==b -> 0.
    {
        const int n = 1;
        std::vector<std::pair<int, int>> edges;
        std::vector<Op> ops{{3, 0, 0}};  // -> 0
        CPF_CHECK(run<BeardHLD>(n, edges, ops) == (std::vector<int>{0}));
        expect_agree(n, edges, ops);
    }

    // a==b on a bigger tree -- always 0, regardless of any edge's color.
    {
        const int n = 5;
        auto edges = path_edges(n);
        std::vector<Op> ops{{2, 0, 0}, {2, 2, 0}, {3, 2, 2}, {3, 0, 0}};
        CPF_CHECK(run<BeardHLD>(n, edges, ops) == (std::vector<int>{0, 0}));
        expect_agree(n, edges, ops);
    }

    // single edge: query it black (1), paint white (-1), paint back black (1).
    {
        const int n = 2;
        std::vector<std::pair<int, int>> edges{{0, 1}};
        std::vector<Op> ops{
            {3, 0, 1},  // black -> 1
            {2, 0, 0},  // paint white
            {3, 0, 1},  // -> -1
            {1, 0, 0},  // paint black again
            {3, 0, 1},  // -> 1
            {3, 1, 0},  // reversed endpoints, same answer -> 1
        };
        CPF_CHECK(run<BeardHLD>(n, edges, ops) ==
                  (std::vector<int>{1, -1, 1, 1}));
        expect_agree(n, edges, ops);
    }

    // a path -- the deepest HLD chains. whiten a middle edge and watch any path
    // that crosses it flip to -1 while paths that avoid it stay black.
    {
        const int n = 7;  // 0-1-2-3-4-5-6, edges 0..5.
        auto edges = path_edges(n);
        std::vector<Op> ops{
            {3, 0, 6},  // whole path black -> 6
            {2, 3, 0},  // paint edge index 3 (joins 3-4) white
            {3, 0, 6},  // crosses it -> -1
            {3, 0, 3},  // stops at vertex 3, avoids it -> 3
            {3, 4, 6},  // lives past it -> 2
            {1, 3, 0},  // paint it back black
            {3, 0, 6},  // -> 6 again
        };
        CPF_CHECK(run<BeardHLD>(n, edges, ops) ==
                  (std::vector<int>{6, -1, 3, 2, 6}));
        expect_agree(n, edges, ops);
    }

    // a star -- every leaf-to-leaf path is exactly two edges through the center.
    // whiten one spoke and only paths using it break.
    {
        const int n = 6;  // center 0, leaves 1..5, edge i-1 = spoke to leaf i.
        auto edges = star_edges(n);
        std::vector<Op> ops{
            {3, 1, 2},  // leaf 1 -> leaf 2, two black spokes -> 2
            {3, 0, 3},  // center -> leaf 3, one spoke -> 1
            {2, 0, 0},  // paint spoke to leaf 1 white (edge index 0)
            {3, 1, 2},  // uses that spoke -> -1
            {3, 2, 3},  // avoids it -> 2
            {3, 0, 1},  // center -> leaf 1 straight over it -> -1
        };
        CPF_CHECK(run<BeardHLD>(n, edges, ops) ==
                  (std::vector<int>{2, 1, -1, 2, -1}));
        expect_agree(n, edges, ops);
    }

    // toggle discipline: white then black on the same edge must restore it, and
    // repeated paints of one color are idempotent.
    {
        const int n = 4;
        auto edges = path_edges(n);
        std::vector<Op> ops{
            {2, 1, 0}, {2, 1, 0},  // whiten edge 1 twice
            {3, 0, 3},             // path crosses it -> -1
            {1, 1, 0}, {1, 1, 0},  // blacken it twice
            {3, 0, 3},             // restored -> 3
        };
        CPF_CHECK(run<BeardHLD>(n, edges, ops) == (std::vector<int>{-1, 3}));
        expect_agree(n, edges, ops);
    }

    // ---- randomized differential vs the path-walk oracle ----------------------
    // thousands of independent programs on random trees. small n so the O(n)-per-
    // query oracle stays cheap; the op mix salts paints of both colors so real
    // lengths and -1 both show up, and every program starts from all-black. any
    // disagreement prints the seed that replays it.
    {
        constexpr int kProgs = 8000;
        cpf::Rng rng(0x0018ULL);
        int diffs = 0;
        std::uint64_t first_bad = 0;

        for (int s = 0; s < kProgs; ++s) {
            const int n = static_cast<int>(rng.in_range(1, 40));
            const auto edges = random_tree(n, rng);
            const int ne = n - 1;

            const int m = static_cast<int>(rng.in_range(1, 60));
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            // skew the paint/query balance per program so some are paint-heavy
            // (more white edges, more -1s) and some are query-heavy.
            const int query_pct = static_cast<int>(rng.in_range(30, 80));
            for (int q = 0; q < m; ++q) {
                const bool is_query = rng.in_range(0, 99) < query_pct;
                if (is_query || ne == 0) {
                    const int a = static_cast<int>(rng.in_range(0, n - 1));
                    const int b = static_cast<int>(rng.in_range(0, n - 1));
                    ops.push_back({3, a, b});
                } else {
                    // lean slightly toward white so -1 answers are plentiful.
                    const int type = rng.in_range(0, 2) == 0 ? 1 : 2;
                    const int e = static_cast<int>(rng.in_range(0, ne - 1));
                    ops.push_back({type, e, 0});
                }
            }

            if (run<BeardHLD>(n, edges, ops) != run<BeardWalk>(n, edges, ops)) {
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

    // ---- larger random trees, still vs the oracle -----------------------------
    // deeper trees mean more chains per path -- the O(log^2 n) walk and the LCA
    // exclusion get exercised harder. oracle stays O(n) per query, so keep counts
    // modest but push n up.
    {
        constexpr int kProgs = 400;
        cpf::Rng rng(0xBEA5DULL);
        int diffs = 0;
        for (int s = 0; s < kProgs; ++s) {
            const int n = static_cast<int>(rng.in_range(50, 500));
            const auto edges = random_tree(n, rng);
            const int ne = n - 1;
            const int m = 300;
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int q = 0; q < m; ++q) {
                const int roll = static_cast<int>(rng.in_range(0, 4));
                if (roll < 2) {  // ~40% queries.
                    const int a = static_cast<int>(rng.in_range(0, n - 1));
                    const int b = static_cast<int>(rng.in_range(0, n - 1));
                    ops.push_back({3, a, b});
                } else {
                    const int type = rng.in_range(0, 1) == 0 ? 1 : 2;
                    const int e = static_cast<int>(rng.in_range(0, ne - 1));
                    ops.push_back({type, e, 0});
                }
            }
            if (run<BeardHLD>(n, edges, ops) != run<BeardWalk>(n, edges, ops)) {
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf("larger-tree differential: %d programs, %d diffs\n", kProgs,
                    diffs);
    }

    return cpf::report();
}
