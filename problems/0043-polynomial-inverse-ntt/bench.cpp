#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

using p0043::i64;
using p0043::kMod;

// the cost the problem charges: newton doubling to precision n = 1e5, each step a
// pair of NTT convolutions. bench the stated worst case. build the input first so
// the timed region is the inverse, not the rng.
int main() {
    cpf::Rng rng(2026);

    constexpr int kN = 100000;
    std::vector<i64> a(static_cast<std::size_t>(kN));
    for (auto& x : a) x = rng.in_range(0, kMod - 1);
    if (a[0] == 0) a[0] = 1;  // a_0 != 0 -- the invertibility precondition.

    // warm once (cold i-cache, first page faults), then time a handful.
    volatile std::size_t sink = 0;
    {
        const std::vector<i64> b = p0043::inverse(a, kN);
        sink += b.size();
    }

    constexpr int kIters = 20;
    auto t0 = std::chrono::steady_clock::now();
    i64 last = 0;
    for (int it = 0; it < kIters; ++it) {
        const std::vector<i64> b = p0043::inverse(a, kN);
        last += b[static_cast<std::size_t>(kN - 1)];
        sink += b.size();
    }
    auto t1 = std::chrono::steady_clock::now();
    cpf::keep(sink);
    cpf::keep(last);

    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / kIters;
    std::printf("polynomial inverse (NTT + newton)\n");
    std::printf("  n=%d\n", kN);
    std::printf("  %8.2f ms/run   %8.2f Mcoeff/s\n", ms,
                static_cast<double>(kN) / ms / 1e3);
    std::printf("  (sink=%zu, last=%lld -- keeps the optimizer honest)\n",
                static_cast<std::size_t>(sink), static_cast<long long>(last));
    return 0;
}
