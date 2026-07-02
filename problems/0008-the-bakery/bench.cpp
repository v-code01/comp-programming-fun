#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case from the constraints: n = 35000, k = 50. that's O(n k log n) window
// steps -- roughly 35000 * 50 * 16 ~ 2.8e7 amortized adds/removes per solve, no
// n^2 anything. types drawn across the full [1, n] alphabet so distinct-counts
// stay lively and the window can't cheat with a tiny type set. pre-generate the
// input so the timed number is the solve, not the rng feeding it.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 35000;
    constexpr int kK = 50;

    std::vector<int> a(kN);
    for (auto& x : a) x = static_cast<int>(rng.in_range(1, kN));

    // one solve is tens of ms, so keep the iteration count small -- the harness
    // warms once and reports per-op, which is the per-solve wall time.
    cpf::bench("bakery solve n=35000 k=50", 3,
               [&] { return p0008::solve(a, kK); });
    std::printf("(ns/op above is the full solve; divide by 1e6 for ms)\n");
    return 0;
}
