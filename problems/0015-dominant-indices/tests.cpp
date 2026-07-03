#include <cstddef>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0015::dominant_indices;
using p0015::dominant_indices_ref;

namespace {

// a tree: vertex count plus its 1-indexed undirected edge list. small enough for
// the O(n^2) oracle to answer instantly.
struct Tree {
    int n;
    std::vector<std::pair<int, int>> edges;
};

// n vertices in a straight line 1-2-...-n. depth 1e6 in the bench, here tiny.
Tree path(int n) {
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) t.edges.emplace_back(i - 1, i);
    return t;
}

// one center (1) with n-1 leaves.
Tree star(int n) {
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) t.edges.emplace_back(1, i);
    return t;
}

// complete binary tree, parent of i is i/2.
Tree balanced_binary(int n) {
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) t.edges.emplace_back(i / 2, i);
    return t;
}

// random rooted tree: vertex i (>= 2) hangs off a uniform earlier vertex. every
// shape from a bamboo to a star shows up as the parent choices vary.
Tree gen(cpf::Rng& rng) {
    const int n = static_cast<int>(rng.in_range(1, 40));
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) {
        const int p = static_cast<int>(rng.in_range(1, i - 1));
        t.edges.emplace_back(p, i);
    }
    return t;
}

}  // namespace

int main() {
    // the two statement shapes, spelled out.
    CPF_EQ(dominant_indices(4, path(4).edges), (std::vector<int>{0, 0, 0, 0}));
    CPF_EQ(dominant_indices(4, star(4).edges), (std::vector<int>{1, 0, 0, 0}));

    // single vertex: no edges, only depth 0.
    CPF_EQ(dominant_indices(1, {}), (std::vector<int>{0}));

    // two vertices: root sees one child at depth 1, count ties depth 0, so 0.
    CPF_EQ(dominant_indices(2, path(2).edges), (std::vector<int>{0, 0}));

    // longer path -- every subtree is a chain, one vertex per level, answer 0.
    CPF_EQ(dominant_indices(6, path(6).edges),
           (std::vector<int>{0, 0, 0, 0, 0, 0}));

    // bigger star -- root's depth 1 dominates with 5, the rest are leaves.
    CPF_EQ(dominant_indices(6, star(6).edges),
           (std::vector<int>{1, 0, 0, 0, 0, 0}));

    // balanced binary tree of 7. root: depth 0 -> 1, depth 1 -> 2, depth 2 -> 4,
    // so the leaves' level wins at 2. each internal-with-two-leaves node peaks at
    // depth 1 (two children beat itself); leaves are 0.
    //   vertices 1..7, parents: 2,3 -> 1; 4,5 -> 2; 6,7 -> 3.
    CPF_EQ(dominant_indices(7, balanced_binary(7).edges),
           (std::vector<int>{2, 1, 1, 0, 0, 0, 0}));

    // a caterpillar to exercise the short-child merges landing on the same depth
    // across several children.  1-2-3 spine, each spine node also carries a leaf.
    //   edges: 1-2, 2-3, 1-4, 2-5, 3-6.
    {
        std::vector<std::pair<int, int>> e{{1, 2}, {2, 3}, {1, 4}, {2, 5}, {3, 6}};
        // hand-check via the oracle so the expectation itself isn't a guess.
        CPF_EQ(dominant_indices(6, e), dominant_indices_ref(6, e));
    }

    // structured shapes over a range of sizes, fast solver against the oracle.
    for (int n = 1; n <= 64; ++n) {
        CPF_EQ(dominant_indices(n, path(n).edges),
               dominant_indices_ref(n, path(n).edges));
        CPF_EQ(dominant_indices(n, star(n).edges),
               dominant_indices_ref(n, star(n).edges));
        CPF_EQ(dominant_indices(n, balanced_binary(n).edges),
               dominant_indices_ref(n, balanced_binary(n).edges));
    }
    std::printf("structured (path/star/binary) n=1..64 vs oracle: ok\n");

    // the offset arithmetic in the shared pool is the bug-prone spot -- so hammer
    // it. thousands of random small trees, every one diffed against the O(n^2)
    // histogram oracle. any disagreement prints the seed that replays it.
    const bool ok = cpf::differential(
        50000, 0xD0'11'CE, gen,
        [](const Tree& t) { return dominant_indices(t.n, t.edges); },
        [](const Tree& t) { return dominant_indices_ref(t.n, t.edges); });
    CPF_CHECK(ok);
    std::printf("differential: 50000 random trees vs O(n^2) oracle\n");

    return cpf::report();
}
