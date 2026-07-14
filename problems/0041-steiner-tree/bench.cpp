#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the solve is one DP over 2^k * n states, so a ns/op number means nothing. what
// matters is the wall for one worst-legal instance: n = 100, k = 10, complete
// graph (m = 4950). that's 1024 * 100 dp cells, 3^10 * 100 merge steps and 1024
// Dijkstras over 4950 edges -- the two halves of the complexity, both maxed.
//
// also time a sparse n = 100, k = 10 graph. same 3^k merge, far cheaper grow, so
// the gap between the two numbers tells you which term is actually paying.
namespace {

double time_once(int n, const std::vector<p0041::Edge>& edges,
                 const std::vector<int>& terminals, std::int64_t* out) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::int64_t ans = p0041::steiner_tree(n, edges, terminals);
    const auto t1 = std::chrono::steady_clock::now();
    cpf::keep(ans);  // pin it -- otherwise -O2 deletes the whole solve.
    *out = ans;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void run(const char* label, int n, const std::vector<p0041::Edge>& edges,
         const std::vector<int>& terminals) {
    std::int64_t ans = 0;
    double best = 1e30;
    for (int trial = 0; trial < 5; ++trial) {
        const double ms = time_once(n, edges, terminals, &ans);
        if (ms < best) best = ms;
        std::printf("  trial %d: %8.3f ms  answer=%lld\n", trial, ms,
                    static_cast<long long>(ans));
    }
    std::printf("%s  n=%d m=%d k=%d  best of 5: %.3f ms  answer=%lld\n\n", label,
                n, static_cast<int>(edges.size()),
                static_cast<int>(terminals.size()), best,
                static_cast<long long>(ans));
}

}  // namespace

int main() {
    cpf::Rng rng(20260713);
    constexpr int kN = 100;
    constexpr int kK = 10;

    // complete graph at the weight ceiling -- the legal maximum on every axis.
    {
        std::vector<p0041::Edge> edges;
        edges.reserve(static_cast<std::size_t>(kN) * (kN - 1) / 2);
        for (int u = 0; u < kN; ++u) {
            for (int v = u + 1; v < kN; ++v) {
                edges.push_back({u, v, rng.in_range(1, 1000000)});
            }
        }
        std::vector<int> terminals;
        for (int i = 0; i < kK; ++i) terminals.push_back(i * 7 % kN);
        run("dense ", kN, edges, terminals);
    }

    // sparse: a random spanning tree plus n extra edges. same k, tiny m.
    {
        std::vector<p0041::Edge> edges;
        for (int v = 1; v < kN; ++v) {
            edges.push_back({static_cast<int>(rng.in_range(0, v - 1)), v,
                             rng.in_range(1, 1000000)});
        }
        for (int i = 0; i < kN; ++i) {
            const int u = static_cast<int>(rng.in_range(0, kN - 1));
            const int v = static_cast<int>(rng.in_range(0, kN - 1));
            if (u != v) edges.push_back({u, v, rng.in_range(1, 1000000)});
        }
        std::vector<int> terminals;
        for (int i = 0; i < kK; ++i) terminals.push_back(i * 9 % kN);
        run("sparse", kN, edges, terminals);
    }

    return 0;
}
