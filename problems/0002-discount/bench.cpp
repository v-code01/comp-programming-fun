#include <cstddef>
#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: n=3e5 costs, m=n-1 coupons -- one of every rank. one op is a full
// solve: sort the n costs, then answer all m coupons at O(1) each. the sort is
// O(n log n) and swamps the answers, so the number is essentially sort-bound.
// inputs are pre-generated -- the timed loop measures the solve, not the rng.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 300000;

    std::vector<std::int64_t> costs = rng.vec(static_cast<std::size_t>(kN), 1,
                                              1000000000);
    std::vector<int> queries;
    queries.reserve(static_cast<std::size_t>(kN - 1));
    for (int q = 2; q <= kN; ++q) queries.push_back(q);  // m = n-1, every rank.

    cpf::bench("discount full solve n=3e5", 200, [&] {
        auto out = p0002::solve(costs, queries);
        return out[out.size() / 2];  // touch the result so the solve can't vanish.
    });
    return 0;
}
