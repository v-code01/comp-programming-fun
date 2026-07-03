#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst shape for the flow: n = 50, q = 100, ranges left wide so the graph is
// dense -- ~50*50 pos->val edges plus 50*50 val->T convex edges, 50 units pushed.
// one solve is far too heavy to loop under cpf::bench, so time a few by hand and
// report the wall clock plus the network size the MCMF chewed through.
int main() {
    constexpr int kN = 50;
    constexpr int kQ = 100;

    // build a feasible dense instance. a couple of band facts squeeze every
    // position into [10, 40] (31 legal values each -- dense), then pad to q=100
    // with redundant loose facts that don't empty any interval. all feasible.
    cpf::Rng rng(0x0020BEEF);
    std::vector<p0020::Fact> facts;
    facts.push_back({1, 1, kN, 10});  // a_i >= 10 everywhere.
    facts.push_back({2, 1, kN, 40});  // a_i <= 40 everywhere.
    while (static_cast<int>(facts.size()) < kQ) {
        int a = static_cast<int>(rng.in_range(1, kN));
        int b = static_cast<int>(rng.in_range(1, kN));
        int l = a < b ? a : b;
        int r = a < b ? b : a;
        // loose bounds inside [10,40] so nothing goes empty: raise floor no
        // higher than 10, drop ceiling no lower than 40.
        if (rng.in_range(0, 1) == 0) {
            facts.push_back({1, l, r, static_cast<int>(rng.in_range(1, 10))});
        } else {
            facts.push_back({2, l, r, static_cast<int>(rng.in_range(40, kN))});
        }
    }

    constexpr int kRuns = 5;
    double best_ms = 1e300;
    std::int64_t answer = 0;
    for (int run = 0; run < kRuns; ++run) {
        auto t0 = std::chrono::steady_clock::now();
        std::int64_t cost = p0020::solve(kN, facts);
        auto t1 = std::chrono::steady_clock::now();
        cpf::keep(cost);  // pin the result so the call can't be elided.
        answer = cost;
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best_ms) best_ms = ms;
    }

    std::printf("dense worst case  n=%d q=%d\n", kN, kQ);
    std::printf("  nodes = 2n+2 = %d, edges = O(n^2), units = n = %d\n",
                2 * kN + 2, kN);
    std::printf("  min cost = %lld\n", static_cast<long long>(answer));
    std::printf("  best of %d runs: %.3f ms\n", kRuns, best_ms);
    return 0;
}
