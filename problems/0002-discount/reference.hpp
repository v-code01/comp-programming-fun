#pragma once

#include <cstdint>
#include <vector>

namespace p0002 {

// independent oracle -- no sort, no order-statistic shortcut. for group size q
// the freed bar is the largest cost c that can be some q-subset's minimum, which
// means at least q bars are priced >= c. count that directly: for each candidate
// cost, sweep the array and tally how many reach it. the largest candidate that
// clears q is the freed value. O(n^2), slow on purpose, right by construction.
// answer = total - freed value.
class ReferenceSolver {
public:
    std::int64_t solve(const std::vector<std::int64_t>& costs, int q) const {
        std::int64_t total = 0;
        for (std::int64_t c : costs) total += c;

        std::int64_t freed = 0;  // costs are >= 1, so 0 is a floor nothing hits.
        for (std::int64_t cand : costs) {
            int atleast = 0;
            for (std::int64_t c : costs)
                if (c >= cand) ++atleast;
            // cand can be a group's minimum only if q bars reach it. keep the max.
            if (atleast >= q && cand > freed) freed = cand;
        }
        return total - freed;
    }
};

}  // namespace p0002
