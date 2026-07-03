#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

// n vertices in one straight chain -- depth n. this is the shape that overflows a
// recursive DFS, and where a per-k greedy done n times goes O(n^2). the iterative
// solver with sqrt-on-the-answer walks it in O(n sqrt n log n) and never touches
// the call stack. proving that at n=1e5 is half the point of the bench.
std::vector<std::pair<int, int>> path_edges(int n) {
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i) e.emplace_back(i - 1, i);
    return e;
}

// one center, n-1 leaves. f collapses to a couple of distinct values, so the
// sqrt loop runs very few greedy passes -- the easy end of the spectrum.
std::vector<std::pair<int, int>> star_edges(int n) {
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i) e.emplace_back(1, i);
    return e;
}

// random rooted tree: vertex i hangs off a uniform earlier vertex. mixed depths
// and fan-outs, the shape a real test file throws at you.
std::vector<std::pair<int, int>> random_edges(int n, std::uint64_t seed) {
    cpf::Rng rng(seed);
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i)
        e.emplace_back(static_cast<int>(rng.in_range(1, i - 1)), i);
    return e;
}

}  // namespace

// worst-case size from the constraints: n = 1e5. three shapes -- the pathological
// chain, a star, and a random tree. one full solve is milliseconds to tens of
// milliseconds, so iterate a few times and read the per-op number as the per-solve
// wall time.
int main() {
    constexpr int kN = 100'000;

    const std::vector<std::pair<int, int>> chain = path_edges(kN);
    cpf::bench("tree-paths path n=1e5", 5,
               [&] { return p0019::you_are_given_a_tree(kN, chain).size(); });

    const std::vector<std::pair<int, int>> st = star_edges(kN);
    cpf::bench("tree-paths star n=1e5", 5,
               [&] { return p0019::you_are_given_a_tree(kN, st).size(); });

    const std::vector<std::pair<int, int>> rnd = random_edges(kN, 2025);
    cpf::bench("tree-paths random n=1e5", 5,
               [&] { return p0019::you_are_given_a_tree(kN, rnd).size(); });

    std::printf("(ns/op above is the full n-answer solve; divide by 1e6 for ms)\n");
    std::printf("path shape ran to completion -- no recursion, no stack overflow\n");
    return 0;
}
