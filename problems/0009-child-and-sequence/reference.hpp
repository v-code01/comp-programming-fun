#pragma once

#include <cstdint>
#include <vector>

namespace p0009 {

// the oracle. a plain array and the definitions read straight off the statement:
// mod loops the range and reduces each cell, sum loops the range and adds, set
// writes one cell. O(n) per op, no cleverness, nothing to get wrong -- that's the
// whole point. the beats tree has to agree with this on every input.
class Reference {
public:
    void build(const std::vector<std::int64_t>& a) { a_ = a; }

    void point_set(int k, std::int64_t x) { a_[static_cast<std::size_t>(k)] = x; }

    void range_mod(int l, int r, std::int64_t x) {
        for (int i = l; i <= r; ++i) a_[static_cast<std::size_t>(i)] %= x;
    }

    std::int64_t range_sum(int l, int r) {
        std::int64_t s = 0;
        for (int i = l; i <= r; ++i) s += a_[static_cast<std::size_t>(i)];
        return s;
    }

private:
    std::vector<std::int64_t> a_;
};

}  // namespace p0009
