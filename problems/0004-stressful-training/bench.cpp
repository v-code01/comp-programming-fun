#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: n = k = 2e5. one full solve is ~58 binary-search probes, each an
// O(n+k) feasibility pass, so the number here is the whole thing end to end.
// pre-generate a feasible board -- big head starts, modest drains -- so the
// search runs to full depth instead of bailing at the first -1. the timed loop
// re-solves the same board; solve() mutates nothing it's handed.
int main() {
    cpf::Rng rng(1132);
    constexpr std::int64_t kN = 200000;
    constexpr std::int64_t kK = 200000;

    std::vector<std::int64_t> a(static_cast<std::size_t>(kN));
    std::vector<std::int64_t> b(static_cast<std::size_t>(kN));
    for (auto& x : a) x = rng.in_range(100000000000LL, 1000000000000LL);  // ~1e11-1e12.
    for (auto& x : b) x = rng.in_range(1, 1000000);                       // ~1e6.

    cpf::bench("solve n=k=2e5", 40, [&] { return p0004::min_charger_power(a, b, kK); });
    return 0;
}
