#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// n = 1e6 random values across the full [-1e9, 1e9] range -- so the heap fills to
// n and the merge fires often. one full solve is the timed unit; ns/op here is
// per-array wall time. pre-generate the input so the timer measures the sweep,
// not the rng.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 1'000'000;
    std::vector<int> a(kN);
    for (auto& x : a)
        x = static_cast<int>(rng.in_range(-1'000'000'000, 1'000'000'000));

    cpf::bench("slope trick n=1e6", 20, [&] { return p0050::solve(a); });
    return 0;
}
