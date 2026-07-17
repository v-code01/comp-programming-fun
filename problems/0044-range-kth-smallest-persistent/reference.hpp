#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0044 {

// the oracle. no persistence, no segment tree, no cleverness -- copy a[l..r] out,
// sort it, hand back the (k-1)-th element. O(len log len) per query, obviously
// correct: the k-th smallest of a range IS the (k-1)-th slot of that range sorted.
// slow on purpose, small only. it exists to disagree with the fast solution on a
// seed we can replay -- nothing more.
class ReferenceKth {
public:
    // keep a copy of the array; queries slice it live.
    void build(const std::vector<std::int64_t>& a) { a_ = a; }

    // k-th smallest in a[l..r], all 1-based, 1 <= k <= r-l+1.
    std::int64_t query(int l, int r, int k) const {
        std::vector<std::int64_t> slice(a_.begin() + (l - 1), a_.begin() + r);
        std::sort(slice.begin(), slice.end());
        return slice[static_cast<std::size_t>(k - 1)];
    }

private:
    std::vector<std::int64_t> a_;
};

}  // namespace p0044
