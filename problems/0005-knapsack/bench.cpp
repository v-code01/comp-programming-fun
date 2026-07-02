#include <array>
#include <cstdint>

#include "common/bench.hpp"
#include "solution.hpp"

using Cnt = std::array<std::int64_t, 9>;

// the stated worst case is W=1e18, every count 1e16. but the pile there weighs
// 3.6e17 < W, so solve() takes the take-all shortcut -- one pass, no dp. the
// number that actually costs something is the large branch: sum > W, which runs
// the full reachability dp over [0, 6720) plus the 840-wide residue sweep. bench
// both so the report shows the shortcut and the real ceiling side by side.
int main() {
    // shortcut case: total < W, answer is the total.
    Cnt take_all{0,
                 10000000000000000LL, 10000000000000000LL, 10000000000000000LL,
                 10000000000000000LL, 10000000000000000LL, 10000000000000000LL,
                 10000000000000000LL, 10000000000000000LL};
    const std::int64_t Wbig = 1000000000000000000LL;

    // large branch: same huge counts, but W under the pile weight so sum > W and
    // the residue+block machinery has to run in full.
    Cnt heavy = take_all;
    const std::int64_t Wmid = 100000000000000000LL;  // 1e17 < 3.6e17.

    cpf::bench("solve take-all (sum<=W)", 2'000'000,
               [&] { return p0005::solve(Wbig, take_all); });
    cpf::bench("solve large branch (dp)", 200'000,
               [&] { return p0005::solve(Wmid, heavy); });
    return 0;
}
