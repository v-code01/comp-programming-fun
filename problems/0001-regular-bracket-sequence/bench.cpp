#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the check is O(1) -- three comparisons, no allocation. so the number here is
// throughput: how many independent queries a second. pre-generate the inputs so
// the timed loop measures the check, not the rng feeding it.
int main() {
    cpf::Rng rng(2025);
    constexpr std::size_t kN = 1u << 16;  // power of two -- index with a mask.
    std::vector<std::array<std::int64_t, 4>> in(kN);
    for (auto& q : in) {
        q = {rng.in_range(0, 1000000000), rng.in_range(0, 1000000000),
             rng.in_range(0, 1000000000), rng.in_range(0, 1000000000)};
    }

    std::size_t i = 0;
    cpf::bench("regular O(1)", 5'000'000, [&] {
        const auto& q = in[i++ & (kN - 1)];
        return p0001::regular(q[0], q[1], q[2], q[3]);
    });
    return 0;
}
