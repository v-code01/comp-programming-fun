#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the build is n-1 max-flows, so a per-op figure means nothing -- what matters is
// the wall to construct the tree at the worst legal size: n = 150, complete graph
// (m = n(n-1)/2 = 11175 edges) with 1e9-scale caps. that's 149 Dinic runs on the
// densest possible input. time the whole construction the way a judge would, then
// show the query side is a table read.
int main() {
    cpf::Rng rng(20250709);
    constexpr int kN = 150;

    std::vector<std::tuple<int, int, std::int64_t>> edges;
    edges.reserve(static_cast<std::size_t>(kN) * (kN - 1) / 2);
    for (int a = 0; a < kN; ++a)
        for (int b = a + 1; b < kN; ++b)
            edges.emplace_back(a, b, rng.in_range(1, 1'000'000'000LL));

    // one honest timed build, repeated to show it's stable and not a fault fluke.
    double best_ms = 1e30;
    std::int64_t checksum = 0;
    for (int trial = 0; trial < 5; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        p0038::GomoryHu tree(kN, edges);
        const auto t1 = std::chrono::steady_clock::now();

        // sum a diagonal of cuts so -O2 can't delete the build.
        checksum = 0;
        for (int i = 0; i + 1 < kN; ++i) checksum += tree.min_cut(i, i + 1);
        cpf::keep(checksum);

        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
        std::printf("trial %d: n=%d m=%d  build %.3f ms  checksum=%lld\n", trial,
                    kN, static_cast<int>(edges.size()), ms,
                    static_cast<long long>(checksum));
    }
    std::printf("best of 5: %.3f ms build  (n=%d, complete graph)\n", best_ms,
                kN);

    // query throughput: after the tree exists, every answer is one array load.
    p0038::GomoryHu tree(kN, edges);
    constexpr std::size_t kMask = (1u << 14) - 1;
    std::vector<std::pair<int, int>> qs(kMask + 1);
    for (auto& q : qs) {
        int u = static_cast<int>(rng.in_range(0, kN - 1));
        int v = static_cast<int>(rng.in_range(0, kN - 1));
        if (u == v) v = (v + 1) % kN;
        q = {u, v};
    }
    std::size_t i = 0;
    cpf::bench("gomory-hu query", 5'000'000, [&] {
        const auto& q = qs[i++ & kMask];
        return tree.min_cut(q.first, q.second);
    });
    return 0;
}
