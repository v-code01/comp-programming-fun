#include <chrono>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the solve is a full max-flow, so a per-op nanosecond figure means nothing --
// what matters is the wall for one worst-legal instance: n = 5000 projects,
// m = 20000 prerequisites. build one random instance at those bounds and time a
// single max_profit call, the way a judge would. pre-generate so the timer
// measures Dinic, not the rng.
int main() {
    cpf::Rng rng(20250703);
    constexpr int kN = 5000;
    constexpr int kM = 20000;

    std::vector<std::int64_t> profits(kN);
    for (auto& p : profits) p = rng.in_range(-1000000000LL, 1000000000LL);

    std::vector<std::pair<int, int>> prereqs;
    prereqs.reserve(kM);
    for (int e = 0; e < kM; ++e) {
        const int i = static_cast<int>(rng.in_range(0, kN - 1));
        const int j = static_cast<int>(rng.in_range(0, kN - 1));
        prereqs.emplace_back(i, j);
    }

    // one honest timed run, plus a couple more to show it's stable and not a
    // page-fault fluke. keep() pins the result so -O2 can't delete the solve.
    std::int64_t answer = 0;
    double best_ms = 1e30;
    for (int trial = 0; trial < 5; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        answer = p0021::max_profit(profits, prereqs);
        const auto t1 = std::chrono::steady_clock::now();
        cpf::keep(answer);
        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
        std::printf("trial %d: n=%d m=%d  %.3f ms  answer=%lld\n", trial, kN, kM,
                    ms, static_cast<long long>(answer));
    }
    std::printf("best of 5: %.3f ms  (n=%d, m=%d)\n", best_ms, kN, kM);
    return 0;
}
