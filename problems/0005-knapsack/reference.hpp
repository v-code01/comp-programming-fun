#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace p0005 {

// the dumb one. a boolean reachability array over [0, W], one item at a time,
// copy by copy. no cleverness to be wrong in -- if dp[v] is set, some subset
// weighs exactly v, and the answer is the last set cell. only sane for small W
// and small counts, which is exactly where the differential test lives.
//
// more than W/i copies of weight i can't appear (one overshoots W), so cap there
// and the loop stays finite even when a count is nominally huge.
inline std::int64_t reference_small(std::int64_t W,
                                    const std::array<std::int64_t, 9>& cnt) {
    std::vector<char> dp(static_cast<std::size_t>(W) + 1, 0);
    dp[0] = 1;
    for (int i = 1; i <= 8; ++i) {
        std::int64_t copies = cnt[i];
        if (copies > W / i) copies = W / i;
        for (std::int64_t c = 0; c < copies; ++c)
            for (std::int64_t v = W; v >= i; --v)
                if (dp[v - i]) dp[v] = 1;
    }
    for (std::int64_t v = W; v >= 0; --v)
        if (dp[v]) return v;
    return 0;
}

// the medium one, and a second independent method on purpose -- the fast solution
// switches strategy above W = 13440, and the small oracle can't reach that far, so
// this carries the large-branch check. 0/1 knapsack over a word-packed bitset:
// split each count into powers of two (1, 2, 4, … items fused into one 0/1 item),
// then OR the shifted bitset in. binary splitting is the textbook exact reduction
// of bounded counts to 0/1 -- obviously correct, just faster than copy-by-copy.
// good up to a few million capacity, which covers the large branch at medium W.
inline std::int64_t reference_medium(std::int64_t W,
                                     const std::array<std::int64_t, 9>& cnt) {
    const std::int64_t words = W / 64 + 1;
    std::vector<std::uint64_t> bits(static_cast<std::size_t>(words), 0);
    bits[0] = 1;  // weight 0 reachable.

    auto add_item = [&](std::int64_t w) {  // OR in bits << w, in place, high->low.
        if (w > W) return;
        const std::int64_t word_sh = w / 64;
        const int bit_sh = static_cast<int>(w % 64);
        for (std::int64_t i = words - 1; i >= 0; --i) {
            std::uint64_t acc = 0;
            std::int64_t src = i - word_sh;
            if (src >= 0) {
                acc = bits[static_cast<std::size_t>(src)] << bit_sh;
                if (bit_sh && src - 1 >= 0)
                    acc |= bits[static_cast<std::size_t>(src - 1)] >> (64 - bit_sh);
            }
            bits[static_cast<std::size_t>(i)] |= acc;
        }
    };

    for (int i = 1; i <= 8; ++i) {
        std::int64_t copies = cnt[i];
        if (copies > W / i) copies = W / i;  // extras overshoot W anyway.
        std::int64_t chunk = 1;
        while (copies > 0) {
            std::int64_t take = chunk < copies ? chunk : copies;
            add_item(take * i);  // one 0/1 item standing in for `take` copies.
            copies -= take;
            chunk <<= 1;
        }
    }

    for (std::int64_t v = W; v >= 0; --v) {
        if (bits[static_cast<std::size_t>(v / 64)] >>
                static_cast<int>(v % 64) & 1u)
            return v;
    }
    return 0;
}

}  // namespace p0005
