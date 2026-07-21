#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case is n=1000, m=1e5. build one dense reachable digraph up front so the
// timed loop measures the solve, not the generation. reachability is guaranteed
// by seeding a root->everyone star with expensive edges -- the random edges then
// fill in the cheap arborescence the algorithm actually chases.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 1000;
    constexpr int kM = 100000;

    std::vector<p0047::Edge> edges;
    edges.reserve(static_cast<std::size_t>(kM));

    // reachability floor: every vertex is one pricey hop from the root.
    for (int v = 2; v <= kN; ++v)
        edges.push_back(p0047::Edge{1, v, 1000000000LL});

    // fill the rest with random directed edges, cheaper so they win the picks.
    while (static_cast<int>(edges.size()) < kM) {
        const int u = static_cast<int>(rng.in_range(1, kN));
        const int v = static_cast<int>(rng.in_range(1, kN));
        const std::int64_t w = rng.in_range(1, 1000000000LL);
        edges.push_back(p0047::Edge{u, v, w});
    }

    std::printf("directed-mst  n=%d  m=%d\n", kN, kM);
    cpf::bench("chu-liu/edmonds", 50, [&] {
        return p0047::min_arborescence(kN, edges, 1);
    });
    return 0;
}
