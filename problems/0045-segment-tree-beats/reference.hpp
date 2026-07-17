#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace p0045 {

// the oracle. a plain array and the four ops read straight off the statement:
// chmin loops the range taking a min, add loops the range adding, max/sum loop
// and fold. O(n) per op, no tree, no lazy, nothing to get wrong -- that is the
// whole point. the beats tree has to agree with this on every input, and when a
// single op-sequence disagrees the tree is what's broken, never this.
class Reference {
public:
    void build(const std::vector<std::int64_t>& a) { a_ = a; }

    void chmin(int l, int r, std::int64_t x) {
        for (int i = l; i <= r; ++i)
            a_[static_cast<std::size_t>(i)] =
                std::min(a_[static_cast<std::size_t>(i)], x);
    }

    void add(int l, int r, std::int64_t x) {
        for (int i = l; i <= r; ++i) a_[static_cast<std::size_t>(i)] += x;
    }

    std::int64_t query_max(int l, int r) {
        std::int64_t m = std::numeric_limits<std::int64_t>::min();
        for (int i = l; i <= r; ++i)
            m = std::max(m, a_[static_cast<std::size_t>(i)]);
        return m;
    }

    std::int64_t query_sum(int l, int r) {
        std::int64_t s = 0;
        for (int i = l; i <= r; ++i) s += a_[static_cast<std::size_t>(i)];
        return s;
    }

private:
    std::vector<std::int64_t> a_;
};

}  // namespace p0045
