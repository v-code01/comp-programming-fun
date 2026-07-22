#include <cstddef>
#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

using p0054::P;

// n=2000 random points on a wide integer grid. concyclic/collinear draws are
// measure-zero at this coord range, so a raw random cloud is general-position
// with overwhelming probability -- no rejection loop needed to time the build.
int main() {
    cpf::Rng rng(20240722);
    const int n = 2000;
    std::vector<P> pts(static_cast<std::size_t>(n));
    for (auto& p : pts) p = {rng.in_range(-10000, 10000), rng.in_range(-10000, 10000)};

    cpf::bench("delaunay n=2000 (Bowyer-Watson O(n^2))", 20,
               [&] { return p0054::delaunay_edge_sq_sum(pts); });
    return 0;
}
