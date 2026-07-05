#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// n = 2e5. the honest worst case for the aliens tree matching: build the tree
// once, then run ~log2(n*maxW) ~ 48 penalized O(n) sweeps inside the toll binary
// search. K = n/4 sits well below the max matching, so the search never collapses
// early. shapes that stress the fold differently:
//   * bamboo  -- a path, deepest possible fold chain, ~n/2 matching.
//   * caterpillar -- a spine with legs, wide fan-in at each spine vertex.
//   * random  -- bushy, mixed depth.
// solve() is one heavy call, not a hot micro-loop, so time it directly with a
// steady clock and print milliseconds.
namespace {

using p0033::Edge;

void time_solve(const char* label, int K, int n, const std::vector<Edge>& e) {
    // warm once -- first touch pays page faults and cold cache.
    p0033::TreeMatch tm(n, e);
    volatile p0033::i64 sink = tm.solve(K);
    (void)sink;

    auto t0 = std::chrono::steady_clock::now();
    const p0033::i64 ans = tm.solve(K);
    auto t1 = std::chrono::steady_clock::now();

    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    cpf::keep(ans);
    std::printf("%-34s n=%d K=%d  ans=%lld  %8.2f ms\n", label, n, K,
                static_cast<long long>(ans), ms);
}

std::vector<Edge> bamboo(int n, cpf::Rng& r) {
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i)
        e.push_back({i, i - 1, r.in_range(1, 1000000000)});
    return e;
}

std::vector<Edge> caterpillar(int n, cpf::Rng& r) {
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    int spine = 1;
    for (int i = 2; i <= n; ++i) {
        const int p = spine;
        if (i % 2 == 0) spine = i;  // even vertices extend the spine.
        e.push_back({i, p, r.in_range(1, 1000000000)});
    }
    return e;
}

std::vector<Edge> random_tree(int n, cpf::Rng& r) {
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 2; i <= n; ++i)
        e.push_back({i, static_cast<int>(r.in_range(1, i - 1)),
                     r.in_range(1, 1000000000)});
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 200000;
    cpf::Rng rng(2033);

    time_solve("bamboo, K=n/4 (binary search)", kN / 4, kN, bamboo(kN, rng));
    time_solve("caterpillar, K=n/4", kN / 4, kN, caterpillar(kN, rng));
    time_solve("random tree, K=n/4", kN / 4, kN, random_tree(kN, rng));

    // a smaller K forces the toll even more negative before the count drops to K,
    // so the search visits the wide negative range -- the deepest part of the log.
    time_solve("bamboo, K=n/8 (deep toll)", kN / 8, kN, bamboo(kN, rng));

    return 0;
}
