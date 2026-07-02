#include <chrono>
#include <cstddef>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "solution.hpp"

// the worst case for O(V*D): min price 1 and every larger price distinct, so the
// shifts b run 1..999 -- D = 999 coins -- and k = 1000 pushes the extra-sum
// ceiling V = k*maxb = 999000. that's ~1e9 inner steps in one call, far too heavy
// to loop under cpf::bench, so time a handful of calls by hand and report the
// wall clock plus the output width the emit loop had to produce.
int main() {
    constexpr int kN = 1000;
    constexpr int kK = 1000;
    std::vector<int> a(static_cast<std::size_t>(kN));
    for (int i = 0; i < kN; ++i) a[static_cast<std::size_t>(i)] = i + 1;  // 1..1000.

    constexpr int kRuns = 3;
    double best_ms = 1e300;
    std::size_t width = 0;
    for (int r = 0; r < kRuns; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<int> out = p0010::solve(kK, a);
        auto t1 = std::chrono::steady_clock::now();
        cpf::keep(out.size());  // pin the result so the call can't be elided.
        width = out.size();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
    }

    std::printf("worst case n=%d k=%d maxcost=%d\n", kN, kK, kN);
    std::printf("  V = k*maxb = %d, D = 999 coins\n", kK * (kN - 1));
    std::printf("  best of %d runs: %.1f ms, %zu totals emitted\n", kRuns,
                best_ms, width);
    return 0;
}
