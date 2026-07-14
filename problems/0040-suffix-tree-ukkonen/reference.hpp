#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace p0040 {

// the oracle. the two definitions, typed out, and nothing else. no tree, no
// automaton, no suffix array -- if it shared machinery with the solution it
// would share the solution's bugs.
//
// cut out every substring s[i..j]. drop them all in a set: the size is (1).
// tally each one in a map: the longest with count >= 2 is (2). occurrences
// overlap freely -- "ana" starts at 1 and at 3 in "banana" and both count,
// which is exactly what a start-position tally gives you.
//
// O(n^2) substrings, each O(n) to cut and hash -- O(n^3), and obvious at a
// glance. small strings only. it's the ground truth the suffix tree has to
// match on every single case.
inline std::pair<std::int64_t, std::int64_t> reference_stats(
    const std::string& s) {
    const int n = static_cast<int>(s.size());
    std::unordered_set<std::string> seen;
    std::unordered_map<std::string, int> count;

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            std::string sub = s.substr(static_cast<std::size_t>(i),
                                       static_cast<std::size_t>(j - i + 1));
            seen.insert(sub);
            ++count[sub];
        }
    }

    std::int64_t longest = 0;
    for (const auto& [sub, c] : count) {
        if (c >= 2) {
            longest = std::max(longest, static_cast<std::int64_t>(sub.size()));
        }
    }
    return {static_cast<std::int64_t>(seen.size()), longest};
}

}  // namespace p0040
