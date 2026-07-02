#pragma once

#include <cstdint>
#include <vector>

namespace p0007 {

// dumb and obviously correct, for small n only. walk all 2^n index-subsets, skip
// the empty one, AND the chosen elements together, and count the subsets that
// land on 0. no cleverness -- this is the definition, spelled out, so the fast
// SOS solution has something independent to disagree with.
//
// n stays <= ~16 in tests, so 2^n subsets is a few tens of thousands. the count
// there is well under the modulus, but take it mod p anyway so the type and the
// value match what solve() returns.
inline std::int64_t brute(const std::vector<std::uint32_t>& a) {
    const int n = static_cast<int>(a.size());
    std::int64_t count = 0;
    for (std::uint32_t sub = 1; sub < (1u << n); ++sub) {
        std::uint32_t andv = ~0u;  // AND identity -- narrows as elements join.
        for (int i = 0; i < n; ++i) {
            if (sub & (1u << i)) andv &= a[i];
        }
        if (andv == 0) ++count;
    }
    return count % 1000000007;
}

}  // namespace p0007
