#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>

namespace p0030 {

// the oracle. the definition, read out loud -- nothing clever.
// enumerate every substring s[i..j], test each for being a palindrome, drop the
// palindromic ones into a set. the set's size is the count of distinct
// palindromic substrings.
//
// O(n^2) substrings, O(n) to test each, plus the substr copy -- O(n^3) and dead
// obvious. only ever run on tiny strings, where it's the ground truth the eertree
// has to match on every case. small alphabets in the differential make palindromes
// dense, so the set fills and the two implementations get stressed where they'd
// disagree.
inline std::int64_t reference_count(const std::string& s) {
    const int n = static_cast<int>(s.size());
    std::unordered_set<std::string> pals;
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            bool is_pal = true;
            for (int a = i, b = j; a < b; ++a, --b) {
                if (s[static_cast<std::size_t>(a)] !=
                    s[static_cast<std::size_t>(b)]) {
                    is_pal = false;
                    break;
                }
            }
            if (is_pal) {
                pals.insert(s.substr(static_cast<std::size_t>(i),
                                     static_cast<std::size_t>(j - i + 1)));
            }
        }
    }
    return static_cast<std::int64_t>(pals.size());
}

}  // namespace p0030
