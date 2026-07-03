#include <cstddef>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0019::you_are_given_a_tree;
using p0019::you_are_given_a_tree_ref;

namespace {

// a tree: vertex count plus its 1-indexed undirected edge list. small enough for
// the exponential oracle to answer instantly.
struct Tree {
    int n;
    std::vector<std::pair<int, int>> edges;
};

// n vertices in a straight line 1-2-...-n. f(k) = floor(n/k) here -- a chain of n
// packs exactly floor(n/k) disjoint length-k paths.
Tree path(int n) {
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) t.edges.emplace_back(i - 1, i);
    return t;
}

// one center (1) with n-1 leaves. longest path is 3 vertices, and every edge
// shares the center -- so f = [n, 1, 1, 0, 0, ...] for n >= 3.
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
    const int n = static_cast<int>(rng.in_range(1, 12));
    Tree t{n, {}};
    for (int i = 2; i <= n; ++i) {
        const int p = static_cast<int>(rng.in_range(1, i - 1));
        t.edges.emplace_back(p, i);
    }
    return t;
}

}  // namespace

int main() {
    // the statement example, spelled out: a 7-vertex path -> "7 3 2 1 1 1 1".
    CPF_EQ(you_are_given_a_tree(7, path(7).edges),
           (std::vector<int>{7, 3, 2, 1, 1, 1, 1}));

    // single vertex: only k=1 exists, one path of one vertex.
    CPF_EQ(you_are_given_a_tree(1, {}), (std::vector<int>{1}));

    // two vertices: k=1 -> 2 singletons, k=2 -> the one edge.
    CPF_EQ(you_are_given_a_tree(2, path(2).edges), (std::vector<int>{2, 1}));

    // a star of 5: 5 singletons, then one edge, then one leaf-center-leaf, then
    // nothing longer than 3 vertices exists.
    CPF_EQ(you_are_given_a_tree(5, star(5).edges),
           (std::vector<int>{5, 1, 1, 0, 0}));

    // balanced binary tree of 7 -- hand-checked against the oracle so the
    // expectation is not itself a guess.
    CPF_EQ(you_are_given_a_tree(7, balanced_binary(7).edges),
           you_are_given_a_tree_ref(7, balanced_binary(7).edges));

    // path shapes obey floor(n/k) exactly -- check the formula directly.
    for (int n = 1; n <= 30; ++n) {
        std::vector<int> want(static_cast<std::size_t>(n));
        for (int k = 1; k <= n; ++k)
            want[static_cast<std::size_t>(k - 1)] = n / k;
        CPF_EQ(you_are_given_a_tree(n, path(n).edges), want);
    }
    std::printf("path floor(n/k) formula n=1..30: ok\n");

    // structured shapes over a range of sizes, fast solver vs the exponential
    // oracle across all k. n kept small so the 2^n memo stays cheap.
    for (int n = 1; n <= 14; ++n) {
        CPF_EQ(you_are_given_a_tree(n, path(n).edges),
               you_are_given_a_tree_ref(n, path(n).edges));
        CPF_EQ(you_are_given_a_tree(n, star(n).edges),
               you_are_given_a_tree_ref(n, star(n).edges));
        CPF_EQ(you_are_given_a_tree(n, balanced_binary(n).edges),
               you_are_given_a_tree_ref(n, balanced_binary(n).edges));
    }
    std::printf("structured (path/star/binary) n=1..14 vs oracle: ok\n");

    // the greedy's optimality is the subtle claim -- so hammer it. thousands of
    // random small trees, each diffed against the independent brute oracle across
    // every k. any disagreement prints the seed that replays it and stops.
    const bool ok = cpf::differential(
        20000, 0x0019BEEF, gen,
        [](const Tree& t) { return you_are_given_a_tree(t.n, t.edges); },
        [](const Tree& t) { return you_are_given_a_tree_ref(t.n, t.edges); });
    CPF_CHECK(ok);
    std::printf("differential: 20000 random trees vs exponential oracle\n");

    return cpf::report();
}
