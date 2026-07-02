#pragma once

#include <array>
#include <cstdint>
#include <unordered_map>

namespace p0001 {

// dumb and obviously correct, for small counts only. try every next block from
// what's left, track the running balance, and refuse to let it go negative --
// including the dip in the middle of a block. memoize on (remaining counts,
// balance) so the search stays finite. if any ordering reaches the end at
// balance zero, a regular sequence exists.
//
// the four blocks as (first half-step, second half-step):
//   "((" = +1,+1   "()" = +1,-1   ")(" = -1,+1   "))" = -1,-1
class ReferenceSolver {
public:
    int solve(std::int64_t c1, std::int64_t c2, std::int64_t c3,
              std::int64_t c4) {
        memo_.clear();
        return can(static_cast<int>(c1), static_cast<int>(c2),
                   static_cast<int>(c3), static_cast<int>(c4), 0)
                   ? 1
                   : 0;
    }

private:
    static constexpr std::array<int, 4> kFirst{+1, +1, -1, -1};
    static constexpr std::array<int, 4> kSecond{+1, -1, +1, -1};

    bool can(int a, int b, int c, int d, int bal) {
        if ((a | b | c | d) == 0) return bal == 0;  // everything placed.
        std::int64_t key = pack(a, b, c, d, bal);
        auto it = memo_.find(key);
        if (it != memo_.end()) return it->second;

        bool ok = false;
        int cnt[4] = {a, b, c, d};
        for (int t = 0; t < 4 && !ok; ++t) {
            if (cnt[t] == 0) continue;
            int mid = bal + kFirst[t];
            int end = mid + kSecond[t];
            if (mid < 0 || end < 0) continue;  // a prefix went negative -- dead.
            --cnt[t];
            ok = can(cnt[0], cnt[1], cnt[2], cnt[3], end);
            ++cnt[t];
        }
        memo_.emplace(key, ok);
        return ok;
    }

    // counts stay under 64 in tests and balance under 256, so one int64 holds
    // the whole state -- keeps the memo a plain hash lookup.
    static std::int64_t pack(int a, int b, int c, int d, int bal) {
        std::int64_t k = a;
        k = k * 64 + b;
        k = k * 64 + c;
        k = k * 64 + d;
        k = k * 256 + (bal + 128);
        return k;
    }

    std::unordered_map<std::int64_t, bool> memo_;
};

}  // namespace p0001
