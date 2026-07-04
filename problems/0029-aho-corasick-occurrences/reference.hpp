#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace p0029 {

// dumb and obviously correct, small inputs only. no automaton, no fail links,
// nothing shared with the solution -- that independence is the point, so a bug
// in one can't hide in the other. for each pattern, slide it across T and count
// every start where all its characters line up; overlaps count because the start
// advances by one, not by the pattern length. sum across all patterns, so a
// repeated pattern is counted once per listing. O(|T| * sum|P|).
inline std::int64_t reference_total(const std::string& text,
                                    const std::vector<std::string>& patterns) {
    const std::size_t n = text.size();
    std::int64_t total = 0;
    for (const std::string& p : patterns) {
        const std::size_t m = p.size();
        if (m == 0 || m > n) continue;  // empty or too long -- zero windows.
        // every start [0 .. n-m] where p matches character for character.
        for (std::size_t start = 0; start + m <= n; ++start) {
            bool hit = true;
            for (std::size_t j = 0; j < m; ++j) {
                if (text[start + j] != p[j]) {
                    hit = false;
                    break;
                }
            }
            if (hit) ++total;
        }
    }
    return total;
}

}  // namespace p0029
