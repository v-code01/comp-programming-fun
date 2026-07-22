#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

using p0052::i64;
using p0052::kMod;

// the stated worst case: n = m = 1e5. build the product tree over m points, one
// root reduction of A, then split down. Horner-per-point would be n*m = 1e10
// multiply-adds here -- the tree is what makes this finish. build the input
// outside the timed region so we measure the eval, not the rng.
int main() {
    cpf::Rng rng(2026);

    constexpr int kN = 100000;
    constexpr int kM = 100000;
    std::vector<i64> a(static_cast<std::size_t>(kN));
    for (auto& x : a) x = rng.in_range(0, kMod - 1);
    std::vector<i64> pts(static_cast<std::size_t>(kM));
    for (auto& x : pts) x = rng.in_range(0, kMod - 1);

    // warm once (cold i-cache, first page faults), then time a handful.
    volatile std::size_t sink = 0;
    {
        const std::vector<i64> r = p0052::multipoint_eval(a, pts);
        sink += r.size();
    }

    constexpr int kIters = 5;
    auto t0 = std::chrono::steady_clock::now();
    i64 last = 0;
    for (int it = 0; it < kIters; ++it) {
        const std::vector<i64> r = p0052::multipoint_eval(a, pts);
        last += r[static_cast<std::size_t>(kM - 1)];
        sink += r.size();
    }
    auto t1 = std::chrono::steady_clock::now();
    cpf::keep(sink);
    cpf::keep(last);

    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / kIters;
    std::printf("multipoint evaluation (NTT product tree + remainder tree)\n");
    std::printf("  n=%d m=%d\n", kN, kM);
    std::printf("  %8.2f ms/run   %8.2f Mpoint/s\n", ms,
                static_cast<double>(kM) / ms / 1e3);
    std::printf("  (sink=%zu, last=%lld -- keeps the optimizer honest)\n",
                static_cast<std::size_t>(sink), static_cast<long long>(last));
    return 0;
}
