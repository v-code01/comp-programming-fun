#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the real workload: n = 1e5 points, full coordinate range. the candidate
// reduction plus four sorts is the whole cost, so this measures O(n log n) end
// to end -- build octant edges, kruskal, done. pre-generate the points so the
// timed loop measures the MST, not the rng.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 100000;
    std::vector<p0049::Point> pts(kN);
    for (auto& p : pts) {
        p.x = rng.in_range(-1000000000, 1000000000);
        p.y = rng.in_range(-1000000000, 1000000000);
    }

    cpf::bench("manhattan_mst n=1e5", 30,
               [&] { return p0049::manhattan_mst(pts); });
    return 0;
}
