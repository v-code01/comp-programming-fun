#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: k=20, so 2^20 ~ 1.05e6 subsets per array. the work is the ranked
// zeta/mobius: 2*(k+1) transforms at O(2^k k) plus the O(2^k k^2) rank
// convolution -- ~2^20 * 20^2 on the order of 4e8 core operations. pre-generate
// the input so the timed number is solve() alone.
//
// one solve allocates a few hundred MB and runs in the tens-to-hundreds of ms
// range, so keep the iteration count low; the harness still warms once and
// reports the per-op wall, which is the per-solve time.
int main() {
    cpf::Rng rng(2025);
    constexpr int kK = 20;
    const std::size_t size = static_cast<std::size_t>(1) << kK;

    std::vector<std::uint32_t> a(size), b(size);
    for (auto& x : a)
        x = static_cast<std::uint32_t>(rng.in_range(0, p0056::kMod - 1));
    for (auto& x : b)
        x = static_cast<std::uint32_t>(rng.in_range(0, p0056::kMod - 1));

    // return one entry so the harness's keep() has a scalar to pin; the whole
    // result vector is still built inside solve(), so nothing is optimized away.
    cpf::bench("subset-sum-conv k=20", 3,
               [&] { return p0056::solve(a, b, kK)[0]; });
    std::printf("(ns/op above is one full k=20 solve; divide by 1e6 for ms)\n");
    return 0;
}
