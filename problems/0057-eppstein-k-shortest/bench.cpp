#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

using p0057::Edge;

// the constraint ceiling: n=1e5, m=5e5, k=1e5. a spanning path 1->2->...->n
// makes everything reachable to t=n, then random cross-edges pile on sidetracks
// and cycles so H_out is fat and H_T melds deep. this is where the persistent
// build and the k-pop best-first both do real work.
std::vector<Edge> big(int n, int m, std::uint64_t seed) {
    cpf::Rng rng(seed);
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(m));
    for (int i = 1; i < n; ++i)
        e.push_back({i, i + 1, rng.in_range(1, 1'000'000'000)});
    while (static_cast<int>(e.size()) < m) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        e.push_back({u, v, rng.in_range(1, 1'000'000'000)});
    }
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 100000;
    constexpr int kM = 500000;
    constexpr int kK = 100000;

    const std::vector<Edge> edges = big(kN, kM, 57);

    // one solve is milliseconds; a few iterations read the per-op number straight
    // as the per-solve wall time. s=1, t=n, so the shortest path is long and the
    // k answers thread the whole graph.
    cpf::bench("eppstein n=1e5 m=5e5 k=1e5", 5, [&] {
        return p0057::k_shortest_walks(kN, 1, kN, kK, edges).size();
    });

    // sanity: print the first three lengths so the run proves it produced real
    // answers, not an empty vector the optimizer could fake.
    const auto ans = p0057::k_shortest_walks(kN, 1, kN, kK, edges);
    std::printf("first three lengths: %lld %lld %lld\n",
                static_cast<long long>(ans[0]), static_cast<long long>(ans[1]),
                static_cast<long long>(ans[2]));
    std::printf("(ns/op above is one full solve; /1e6 = ms)\n");
    return 0;
}
