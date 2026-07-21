#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

using Edge = std::tuple<int, int, int>;

// the worst case the constraints allow: n = m = 300. random edges over 300
// vertices, colors spread across 250 of the palette so the partition cap binds
// yet still admits a large forest -- the number of augmentations climbs toward
// r = min(n-1, #colors) and every step pays a full O(m*n) exchange-graph build.
std::vector<Edge> big_graph(int n, int m, int ncolors, std::uint64_t seed) {
    cpf::Rng rng(seed);
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        const int c = static_cast<int>(rng.in_range(1, ncolors));
        e.emplace_back(u, v, c);
    }
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 300, kM = 300;
    const std::vector<Edge> g = big_graph(kN, kM, 250, 2048);

    // one solve is well under a millisecond; iterate enough to read a stable
    // per-solve wall time straight off the ns/op.
    cpf::bench("colorful-forest n=m=300", 2000, [&] {
        return p0048::max_colorful_forest(kN, g);
    });

    std::printf("(ns/op above is one full n=m=300 solve)\n");
    return 0;
}
