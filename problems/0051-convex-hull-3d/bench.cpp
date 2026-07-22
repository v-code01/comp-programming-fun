#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the real workload: n = 2000 points across the full coordinate range. the hull
// is O(n^2) -- every insertion scans O(n) faces -- so this times the whole
// incremental build, ~8e6 __int128 orientation calls end to end. random points
// at 21-bit coords are general-position with overwhelming odds; we don't gate on
// it here because timing doesn't care, only the sign path does. pre-generate the
// points so the timed loop measures the hull, not the rng.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 2000;
    std::vector<p0051::P> pts(static_cast<std::size_t>(kN));
    for (auto& p : pts) {
        p.x = rng.in_range(-1000000, 1000000);
        p.y = rng.in_range(-1000000, 1000000);
        p.z = rng.in_range(-1000000, 1000000);
    }

    cpf::bench("convex_hull_faces n=2000", 50,
               [&] { return p0051::convex_hull_faces(pts); });
    return 0;
}
