#include <cstddef>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

// n vertices in one straight chain -- depth n. this is the shape that overflows a
// recursive DFS, and the shape where a naive per-vertex merge goes O(n^2). the
// iterative long-path solver walks it in O(n) with a bounded stack. proving that
// here is half the point of the bench.
std::vector<std::pair<int, int>> path_edges(int n) {
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i) e.emplace_back(i - 1, i);
    return e;
}

// a random rooted tree: vertex i hangs off a uniform earlier vertex. mixed depths
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

// worst-case size from the constraints: n = 1e6. two shapes -- the pathological
// chain and a random tree. one solve is milliseconds, so iterate a few times and
// read the per-op number as the per-solve wall time.
int main() {
    constexpr int kN = 1'000'000;

    const std::vector<std::pair<int, int>> chain = path_edges(kN);
    cpf::bench("dominant path n=1e6", 3,
               [&] { return p0015::dominant_indices(kN, chain).size(); });

    const std::vector<std::pair<int, int>> rnd = random_edges(kN, 2025);
    cpf::bench("dominant random n=1e6", 3,
               [&] { return p0015::dominant_indices(kN, rnd).size(); });

    std::printf("(ns/op above is the full solve; divide by 1e6 for ms)\n");
    std::printf("path shape ran to completion -- no recursion, no stack overflow\n");
    return 0;
}
