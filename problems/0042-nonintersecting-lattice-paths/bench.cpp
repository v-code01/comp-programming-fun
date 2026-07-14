#include <chrono>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "solution.hpp"

namespace {

// the worst case the constraints allow: k = 100, coordinates spread across the
// full 0..2e5 range, and every end reachable from every start. that last part is
// what makes it worst -- a dense M with no structural zeros means the elimination
// never bails early and every one of the k^3 inner steps runs.
//
// starts walk south-east from (0, 1e5) to (9900, 90100). ends walk south-east from
// (1e5, 2e5) down to (109900, 190100). every b_j is north-east of every a_i, so all
// 10000 entries of M are nonzero binomials.
void worst_case(std::vector<p0042::Pt>& a, std::vector<p0042::Pt>& b) {
    constexpr int k = 100;
    a.resize(k);
    b.resize(k);
    for (int i = 0; i < k; ++i) {
        a[static_cast<std::size_t>(i)] = p0042::Pt{i * 100, 100000 - i * 100};
        b[static_cast<std::size_t>(i)] = p0042::Pt{100000 + i * 100, 200000 - i * 100};
    }
}

}  // namespace

int main() {
    // the factorial tables are O(maxCoord) and built once. time that separately --
    // folding a 4e5 table build into a k^3 determinant would slander the
    // determinant.
    const auto t0 = std::chrono::steady_clock::now();
    p0042::prepare();
    const auto t1 = std::chrono::steady_clock::now();
    std::printf("factorial tables to 4e5 (one time): %.2f ms\n",
                std::chrono::duration<double, std::milli>(t1 - t0).count());

    std::vector<p0042::Pt> a, b;
    worst_case(a, b);

    cpf::bench("LGV k=100 dense, coords 2e5", 200, [&] { return p0042::solve(a, b); });
    std::printf("(ns/op above is one full solve: 10000 binomials + a 100x100 determinant.)\n");
    std::printf("k=100 answer mod 1e9+7 = %lld\n",
                static_cast<long long>(p0042::solve(a, b)));
    return 0;
}
