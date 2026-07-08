#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

using Edges = std::vector<std::pair<int, int>>;

// a 2e5-deep chain 1->2->...->n. this is the case a recursive DFS dies on and
// the one where the link-eval ancestor chain grows longest, so path compression
// earns its keep. n-1 edges, so pad to the m budget with random cross-edges that
// don't shorten anyone's dominator chain much.
Edges chain_plus_cross(int n, int m, std::uint64_t seed) {
    cpf::Rng rng(seed);
    Edges e;
    e.reserve(static_cast<std::size_t>(m));
    for (int i = 1; i < n; ++i) e.emplace_back(i, i + 1);
    while (static_cast<int>(e.size()) < m) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        e.emplace_back(u, v);
    }
    return e;
}

// a wide random digraph: every vertex reachable, heavy merging, no structure the
// solver can shortcut. the general worst case at the constraint ceiling.
Edges random_dense(int n, int m, std::uint64_t seed) {
    cpf::Rng rng(seed);
    Edges e;
    e.reserve(static_cast<std::size_t>(m));
    // a spanning path first so the whole graph is reachable from vertex 1.
    for (int i = 1; i < n; ++i) e.emplace_back(i, i + 1);
    while (static_cast<int>(e.size()) < m) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        e.emplace_back(u, v);
    }
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 200000;
    constexpr int kM = 300000;

    const Edges chain = chain_plus_cross(kN, kM, 36);
    const Edges dense = random_dense(kN, kM, 3600);

    // one solve is milliseconds; a few iterations read the per-op number straight
    // as the per-solve wall time.
    cpf::bench("dominator chain n=2e5 m=3e5", 10, [&] {
        return p0036::dominator_tree(kN, chain, 1).size();
    });
    cpf::bench("dominator dense n=2e5 m=3e5", 10, [&] {
        return p0036::dominator_tree(kN, dense, 1).size();
    });

    std::printf("(ns/op above is one full n=2e5 solve; /1e6 = ms)\n");
    return 0;
}
