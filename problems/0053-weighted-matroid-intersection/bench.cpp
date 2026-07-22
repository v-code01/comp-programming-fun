#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

using ll = std::int64_t;
using Edge = std::tuple<int, int, int, ll>;

// the worst case the constraints allow: n = m = 200. random edges over 200
// vertices, colors spread across 150 of the palette so the cap binds yet still
// admits a large forest -- the augmentation count climbs toward r = min(n-1,
// #colors), and every step pays a full exchange-graph build plus Bellman-Ford.
// mixed-sign weights so the concave stop can actually fire mid-run.
std::vector<Edge> big_graph(int n, int m, int ncolors, std::uint64_t seed) {
    cpf::Rng rng(seed);
    std::vector<Edge> e;
    e.reserve(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        const int c = static_cast<int>(rng.in_range(1, ncolors));
        const ll w = rng.in_range(-1000000, 1000000);
        e.emplace_back(u, v, c, w);
    }
    return e;
}

}  // namespace

int main() {
    constexpr int kN = 200, kM = 200;
    const std::vector<Edge> g = big_graph(kN, kM, 150, 2048);

    cpf::bench("wmi-colorful-forest n=m=200", 1000, [&] {
        return p0053::max_weight_colorful_forest(kN, g);
    });

    std::printf("(ns/op above is one full n=m=200 solve)\n");
    return 0;
}
