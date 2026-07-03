#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0017::Reference;
using p0017::Tree;

// one query in the statement's alphabet.
//   type 1: repaint(v, c)   -- paint v's subtree color c
//   type 2: distinct(v)     -- count colors in v's subtree (c unused)
struct Op {
    int type;
    int v;
    int c;
};

// replay a program on anything with the Tree API, collecting the type-2 answers
// in order. templated so the segment tree and the oracle run the identical path.
template <typename DS>
std::vector<int> run(DS& ds, int n, const std::vector<int>& color,
                     const std::vector<std::pair<int, int>>& edges,
                     const std::vector<Op>& ops) {
    ds.build(n, color, edges);
    std::vector<int> out;
    for (const Op& o : ops) {
        if (o.type == 1)
            ds.repaint(o.v, o.c);
        else
            out.push_back(ds.distinct(o.v));
    }
    return out;
}

void expect_agree(int n, const std::vector<int>& color,
                  const std::vector<std::pair<int, int>>& edges,
                  const std::vector<Op>& ops) {
    Tree tree;
    Reference ref;
    CPF_CHECK(run(tree, n, color, edges, ops) ==
              run(ref, n, color, edges, ops));
}

// a random rooted tree: vertex i (>= 2) hangs off a uniformly random earlier
// vertex. that alone reaches every shape -- a fresh random parent for each node
// spans stars, paths, and everything between.
std::vector<std::pair<int, int>> random_tree(cpf::Rng& rng, int n) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 2; i <= n; ++i) {
        int p = static_cast<int>(rng.in_range(1, i - 1));
        edges.emplace_back(p, i);
    }
    return edges;
}

}  // namespace

int main() {
    // ---- the brute-verified example (Codeforces 620E sample) ---------------
    // 7 vertices, all color 1. tree: 1-{2,3,4}, 3-{5,6,7}. traced by hand and
    // confirmed against the oracle below.
    {
        int n = 7;
        std::vector<int> color{0, 1, 1, 1, 1, 1, 1, 1};  // 1-indexed.
        std::vector<std::pair<int, int>> edges{
            {1, 2}, {1, 3}, {1, 4}, {3, 5}, {3, 6}, {3, 7}};
        std::vector<Op> ops{
            {1, 3, 2},  // subtree(3)={3,5,6,7} -> color 2
            {2, 1, 0},  // colors in whole tree: {1,2}          -> 2
            {1, 4, 3},  // subtree(4)={4} -> color 3
            {2, 1, 0},  // {1,2,3}                               -> 3
            {1, 2, 5},  // subtree(2)={2} -> color 5
            {2, 1, 0},  // {1,2,3,5}                             -> 4
            {1, 6, 4},  // subtree(6)={6} -> color 4
            {2, 1, 0},  // {1,2,3,4,5}                           -> 5
            {2, 2, 0},  // subtree(2)={2}, color 5               -> 1
            {2, 3, 0},  // subtree(3)={3,5,6,7} colors {2,4}     -> 2
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        std::vector<int> want{2, 3, 4, 5, 1, 2};
        CPF_CHECK(got == want);
        expect_agree(n, color, edges, ops);
    }

    // ---- enumerated edges --------------------------------------------------

    // single node: its subtree is itself. repaint then query is trivial but the
    // segment tree of size 1 is a real edge (cap = 4, one leaf).
    {
        int n = 1;
        std::vector<int> color{0, 7};
        std::vector<std::pair<int, int>> edges;
        std::vector<Op> ops{{2, 1, 0}, {1, 1, 42}, {2, 1, 0}};
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{1, 1}));  // one color, then one color.
        expect_agree(n, color, edges, ops);
    }

    // a path 1-2-3-4-5: deepest possible tree for its size, the shape that
    // stresses the iterative tour. subtree(v) is the tail from v down.
    {
        int n = 5;
        std::vector<int> color{0, 1, 2, 3, 4, 5};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 3}, {3, 4}, {4, 5}};
        std::vector<Op> ops{
            {2, 1, 0},   // whole path, 5 colors                -> 5
            {2, 3, 0},   // subtree(3)={3,4,5} colors {3,4,5}   -> 3
            {1, 2, 9},   // subtree(2)={2,3,4,5} -> color 9
            {2, 1, 0},   // {1,9}                               -> 2
            {2, 3, 0},   // {9}                                 -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{5, 3, 2, 1}));
        expect_agree(n, color, edges, ops);
    }

    // a star: root 1 with n-1 leaves. subtree(root) is everything; each leaf is
    // alone. repaint a leaf touches exactly one cell.
    {
        int n = 5;
        std::vector<int> color{0, 1, 2, 2, 2, 2};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}, {1, 4}, {1, 5}};
        std::vector<Op> ops{
            {2, 1, 0},   // {1,2}                               -> 2
            {1, 3, 3},   // leaf 3 -> color 3
            {2, 1, 0},   // {1,2,3}                             -> 3
            {2, 3, 0},   // just leaf 3                         -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{2, 3, 1}));
        expect_agree(n, color, edges, ops);
    }

    // repaint the whole tree then query the root: everything collapses to one
    // color -- the lazy assignment has to overwrite every cell.
    {
        int n = 6;
        std::vector<int> color{0, 1, 2, 3, 4, 5, 6};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3},
                                               {2, 4}, {2, 5}, {3, 6}};
        std::vector<Op> ops{
            {2, 1, 0},   // 6 distinct colors                  -> 6
            {1, 1, 4},   // repaint the whole tree color 4
            {2, 1, 0},   // one color everywhere               -> 1
            {2, 2, 0},   // subtree(2), also one color         -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{6, 1, 1}));
        expect_agree(n, color, edges, ops);
    }

    // nested repaints: an outer subtree assign, then an inner one carves a
    // different color back out of it -- the lazy tag must push before the inner
    // window splits the segment.
    {
        int n = 5;
        std::vector<int> color{0, 1, 1, 1, 1, 1};
        std::vector<std::pair<int, int>> edges{{1, 2}, {2, 3}, {3, 4}, {4, 5}};
        std::vector<Op> ops{
            {1, 2, 5},   // subtree(2)={2,3,4,5} -> color 5
            {1, 4, 6},   // inner subtree(4)={4,5} -> color 6
            {2, 1, 0},   // v1=1, {2,3}=5, {4,5}=6 -> {1,5,6}  -> 3
            {2, 2, 0},   // {5,6}                              -> 2
            {2, 4, 0},   // {6}                                -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{3, 2, 1}));
        expect_agree(n, color, edges, ops);
    }

    // all-same-color to start: every query is 1 until a repaint introduces a
    // second color.
    {
        int n = 7;
        std::vector<int> color{0, 3, 3, 3, 3, 3, 3, 3};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}, {2, 4},
                                               {2, 5}, {3, 6}, {3, 7}};
        std::vector<Op> ops{
            {2, 1, 0},   // all color 3                        -> 1
            {2, 2, 0},   // all color 3                        -> 1
            {1, 5, 8},   // one leaf -> color 8
            {2, 1, 0},   // {3,8}                              -> 2
            {2, 3, 0},   // subtree(3) untouched, all 3        -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{1, 1, 2, 1}));
        expect_agree(n, color, edges, ops);
    }

    // color at the ceiling: 60 is the top allowed color, bit 59 of the mask.
    // exercise the high bit so an off-by-one in 1<<(c-1) would show.
    {
        int n = 3;
        std::vector<int> color{0, 60, 1, 60};
        std::vector<std::pair<int, int>> edges{{1, 2}, {1, 3}};
        std::vector<Op> ops{
            {2, 1, 0},   // {60,1}                             -> 2
            {1, 1, 60},  // all -> 60
            {2, 1, 0},   // {60}                               -> 1
        };
        Tree tree;
        std::vector<int> got = run(tree, n, color, edges, ops);
        CPF_CHECK((got == std::vector<int>{2, 1}));
        expect_agree(n, color, edges, ops);
    }

    // ---- randomized differential vs the oracle -----------------------------
    // thousands of independent programs on random trees. small n so the O(n)-per-
    // op oracle stays cheap; a tiny color alphabet so repaints collide and the
    // distinct-counts actually move; a repaint-heavy op mix so the lazy assign
    // path gets hammered.
    {
        constexpr int kSeqs = 8000;
        cpf::Rng rng(0x0017ULL);
        int diffs = 0;
        std::uint64_t first_bad = 0;

        for (int s = 0; s < kSeqs; ++s) {
            int n = static_cast<int>(rng.in_range(1, 40));
            // small alphabet forces collisions; occasionally widen toward the
            // 60-color ceiling so the high bits get touched too.
            int chi = rng.in_range(0, 3) == 0 ? 60 : 6;

            std::vector<int> color(static_cast<std::size_t>(n) + 1, 0);
            for (int i = 1; i <= n; ++i)
                color[static_cast<std::size_t>(i)] =
                    static_cast<int>(rng.in_range(1, chi));

            std::vector<std::pair<int, int>> edges = random_tree(rng, n);

            int m = static_cast<int>(rng.in_range(1, 60));
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int q = 0; q < m; ++q) {
                int v = static_cast<int>(rng.in_range(1, n));
                if (rng.in_range(0, 2) == 0) {
                    ops.push_back({2, v, 0});  // ~1/3 queries
                } else {
                    int c = static_cast<int>(rng.in_range(1, chi));
                    ops.push_back({1, v, c});  // ~2/3 repaints
                }
            }

            Tree tree;
            Reference ref;
            if (run(tree, n, color, edges, ops) !=
                run(ref, n, color, edges, ops)) {
                if (diffs == 0) first_bad = static_cast<std::uint64_t>(s);
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs)
            std::printf("  first differing sequence index: %llu\n",
                        static_cast<unsigned long long>(first_bad));
        std::printf("differential: %d programs, %d diffs\n", kSeqs, diffs);
    }

    return cpf::report();
}
