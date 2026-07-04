#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// bellman-ford is O(n * m), so a per-op nanosecond figure means nothing -- what
// matters is the wall for one worst-legal instance. n = 300 currencies, m = n*n =
// 90000 offers, the densest graph the bounds allow. an injected negative cycle
// keeps every sweep productive, so the solve runs all n-1 sweeps -- the true
// worst case, no early bail. pre-build so the timer measures the detect, not the
// rng.
int main() {
    cpf::Rng rng(20260704);
    constexpr int kN = 300;
    const int kM = kN * kN;

    std::vector<p0026::Offer> offers;
    offers.reserve(static_cast<std::size_t>(kM));
    for (int u = 1; u <= kN; ++u)
        for (int v = 1; v <= kN; ++v) {
            // log-uniform rate in [~0.05, ~20] -- spread across the -log axis,
            // all inside [1e-6, 1e6].
            const double lg = static_cast<double>(rng.in_range(-300, 300)) / 100.0;
            offers.push_back({u, v, std::exp(lg)});
        }
    // bolt on one loud loop so the answer is YES and no sweep gets to bail early.
    offers.push_back({1, 2, 1.2});
    offers.push_back({2, 3, 1.2});
    offers.push_back({3, 1, 1.2});

    bool answer = false;
    double best_ms = 1e30;
    for (int trial = 0; trial < 5; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        answer = p0026::has_arbitrage(kN, offers);
        const auto t1 = std::chrono::steady_clock::now();
        cpf::keep(answer);
        const double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
        std::printf("trial %d: n=%d m=%d  %.3f ms  arbitrage=%s\n", trial, kN, kM,
                    ms, answer ? "YES" : "NO");
    }
    std::printf("best of 5: %.3f ms  (n=%d, m=%d)\n", best_ms, kN, kM);
    return 0;
}
