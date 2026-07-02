#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: n = 1e6 elements over the full 20-bit domain. the work is the
// superset SOS -- 2^20 * 20 ~ 21M folds -- plus one O(n) pass to bucket and one
// to sum powers. pre-generate the input so the timed number is the solve alone.
//
// one solve is on the order of milliseconds, so keep the iteration count small;
// the harness still warms once and reports per-op, which is the per-solve wall.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 1'000'000;
    constexpr int kBits = 20;
    const std::uint32_t hi = (1u << kBits) - 1;

    std::vector<std::uint32_t> a(kN);
    for (auto& x : a) x = static_cast<std::uint32_t>(rng.in_range(0, hi));

    cpf::bench("jzzhu solve n=1e6 b=20", 20,
               [&] { return p0007::solve(a, kBits); });
    std::printf("(ns/op above is the full solve; divide by 1e6 for ms)\n");
    return 0;
}
