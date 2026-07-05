#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

using ll = std::int64_t;
using Edge = std::tuple<int, int, ll>;

// the worst case the constraints allow: n = 500, the COMPLETE graph. every one
// of the n(n-1)/2 ~ 124,750 pairs carries a weight, so the blossom scan can't
// skip a thing and the O(n^3) core runs flat out. random weights to keep the
// dual adjustments from collapsing to a trivial schedule.
std::vector<Edge> complete_graph(int n, std::uint64_t seed) {
    cpf::Rng rng(seed);
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n - 1) / 2);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            e.emplace_back(i, j, rng.in_range(1, 1000000000LL));
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 500;
    const std::vector<Edge> dense = complete_graph(kN, 2032);

    // one solve is tens to hundreds of ms; a handful of iterations reads the
    // per-op number straight as the per-solve wall time.
    cpf::bench("max-weight-matching K500", 3, [&] {
        return p0032::max_weight_matching(kN, dense);
    });

    std::printf("(ns/op above is one full n=500 complete-graph solve; /1e6 = ms)\n");
    return 0;
}
