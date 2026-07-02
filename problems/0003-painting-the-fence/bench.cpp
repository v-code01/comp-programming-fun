#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: n = q = 5000. the whole solve is the difference-array sweep plus
// the sparse pair-min -- O(n + q log q), no q*q anything. pre-generate the input
// so the timed number is the solve, not the rng feeding it.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 5000;
    constexpr int kQ = 5000;

    std::vector<std::pair<int, int>> ranges(kQ);
    for (auto& p : ranges) {
        int l = static_cast<int>(rng.in_range(1, kN));
        int r = static_cast<int>(rng.in_range(1, kN));
        if (l > r) std::swap(l, r);
        p = {l, r};
    }

    // one solve is tens of ms, so keep the iteration count tiny -- the harness
    // still warms once and reports per-op, which is the per-solve wall time.
    cpf::bench("painting solve n=q=5000", 3,
               [&] { return p0003::solve(kN, kQ, ranges); });
    std::printf("(ns/op above is the full solve; divide by 1e6 for ms)\n");
    return 0;
}
