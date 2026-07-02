#include <cstddef>
#include <cstdint>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case: n = 1e5. a strictly increasing from 1 up toward 1e9, b strictly
// decreasing down to 0 -- random steps so the hull churns (lines get inserted and
// retired) instead of collapsing to one flat envelope. one min_cost is O(n), so
// the number here is a full solve end to end. pre-build the board; min_cost reads
// it and mutates nothing, so re-solving the same one each iteration is honest.
int main() {
    cpf::Rng rng(319);
    constexpr std::int64_t kN = 100000;

    std::vector<std::int64_t> a(static_cast<std::size_t>(kN));
    std::vector<std::int64_t> b(static_cast<std::size_t>(kN));

    // steps in [1, 9000]: n * 9000 <= 9e8 < 1e9, so a[n-1] and b[0] stay under
    // the 1e9 bound while both sequences stay strictly monotone.
    a[0] = 1;
    for (std::size_t i = 1; i < a.size(); ++i) a[i] = a[i - 1] + rng.in_range(1, 9000);
    b[b.size() - 1] = 0;
    for (std::size_t i = b.size() - 1; i-- > 0;) b[i] = b[i + 1] + rng.in_range(1, 9000);

    cpf::bench("min_cost n=1e5", 2000, [&] { return p0011::min_cost(a, b); });
    return 0;
}
