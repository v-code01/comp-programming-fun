#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>

namespace p0012 {

// the oracle. no automaton, no suffix links -- just the definition read out loud.
// a rotation of x is x+x sliced at each of the |x| offsets; collect the distinct
// rotation-strings into a set. then slide a length-|x| window across s and count
// the windows that land in that set. each window is one position, counted once,
// so rotations that coincide never inflate the total.
//
// O(|s| * |x|) for the scan (a substr + hash lookup per window) plus O(|x|^2) to
// build the rotation set. dead simple, obviously right -- exactly what a fast
// solution has to agree with on every small case, especially the ones where
// rotations collide and the dedupe has to earn its keep.
inline std::int64_t reference_count(const std::string& s, const std::string& x) {
    const int L = static_cast<int>(x.size());
    const int n = static_cast<int>(s.size());
    if (L == 0 || L > n) return 0;

    const std::string doubled = x + x;
    std::unordered_set<std::string> rotations;
    rotations.reserve(static_cast<std::size_t>(L) * 2);
    for (int k = 0; k < L; ++k) {
        rotations.insert(doubled.substr(static_cast<std::size_t>(k),
                                        static_cast<std::size_t>(L)));
    }

    std::int64_t ans = 0;
    for (int i = 0; i + L <= n; ++i) {
        if (rotations.count(s.substr(static_cast<std::size_t>(i),
                                     static_cast<std::size_t>(L))) != 0) {
            ++ans;
        }
    }
    return ans;
}

}  // namespace p0012
