#pragma once

#include <cstdint>
#include <vector>

namespace p0007 {

constexpr std::int64_t kMod = 1000000007;

// count the non-empty index-subsets whose bitwise AND is exactly 0, mod 1e9+7.
//
// the honest count is over 2^n subsets -- dead at n = 1e6. so flip the question.
// "AND == 0" is 20 independent "bit b is clear in the AND" conditions, and a bit
// is clear in the AND iff at least one chosen element has it clear. that's a
// union of bad events, one per bit … which is what inclusion-exclusion eats.
//
// fix a set S of bits and force all of them SET in the AND. the AND carries bit b
// iff every chosen element carries bit b, so forcing S set means every element in
// the subset is a supermask of S. let c(S) = # of the a_i that are supermasks of
// S. the non-empty subsets drawn only from those elements number 2^c(S) - 1, and
// every one of them has an AND that is a supermask of S.
//
// by inclusion-exclusion over which bits the AND is forced to carry:
//   answer = sum over S of (-1)^popcount(S) * (2^c(S) - 1).
// the "- 1" drops the empty subset in every term, so the empty group never leaks
// in -- S = 0 contributes 2^n - 1 (all non-empty subsets) and the higher terms
// carve out the ones whose AND kept a bit.
//
// c(S) is a superset sum of the frequency array: c(S) = sum over T superset of S
// of freq(T). the standard SOS transform does all 2^B of them in O(2^B * B):
// for each bit, fold the "bit set" half into the "bit clear" half.
//
// bits = number of low bits in play (20 for the real input). every a_i must be
// < 2^bits. complexity O(2^bits * bits + n), space O(2^bits + n).
inline std::int64_t solve(const std::vector<std::uint32_t>& a, int bits) {
    const std::size_t n = a.size();
    const std::size_t size = static_cast<std::size_t>(1) << bits;

    // freq(v) = how many a_i equal v. counts stay <= n <= 1e6, an int holds them.
    std::vector<std::int32_t> c(size, 0);
    for (std::uint32_t v : a) ++c[v];

    // superset-sum SOS: after bit i is folded, c(S) counts every element that
    // agrees with S on the bits seen so far and is a supermask on bit i. run all
    // bits and c(S) is the full supermask count. c(S) <= n, no overflow.
    for (int i = 0; i < bits; ++i) {
        const std::size_t hi = static_cast<std::size_t>(1) << i;
        for (std::size_t mask = 0; mask < size; ++mask) {
            if ((mask & hi) == 0) c[mask] += c[mask | hi];
        }
    }

    // powers of two mod p, up to n -- c(S) never exceeds n.
    std::vector<std::int64_t> pow2(n + 1);
    pow2[0] = 1;
    for (std::size_t i = 1; i <= n; ++i) pow2[i] = pow2[i - 1] * 2 % kMod;

    // sum the signed terms. popcount parity is the sign; keep it non-negative by
    // adding p before every subtract so the running total stays in [0, p).
    std::int64_t ans = 0;
    for (std::size_t mask = 0; mask < size; ++mask) {
        std::int64_t term = (pow2[c[mask]] - 1 + kMod) % kMod;
        if (__builtin_popcountll(mask) & 1) {
            ans = (ans - term + kMod) % kMod;
        } else {
            ans = (ans + term) % kMod;
        }
    }
    return ans;
}

}  // namespace p0007
