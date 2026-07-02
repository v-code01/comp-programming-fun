#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace p0002 {

// coupon q: buy exactly q bars in one purchase and the cheapest of those q is
// free. you buy all n bars, and you pick which q ride the coupon. to pay the
// least you want the freed bar -- the group's minimum -- as expensive as it can
// be. the min of any q-subset is largest when you take the q priciest bars:
// their minimum is the q-th largest cost in the whole set. so
//   answer(q) = total - (q-th largest cost).
// the other n-q bars pay face value and their sum is fixed -- they never enter
// the choice. one descending sort lays out every order statistic at once, so
// each coupon answers in O(1) after it.
class Discount {
public:
    explicit Discount(std::vector<std::int64_t> costs)
        : desc_(std::move(costs)), total_(0) {
        for (std::int64_t c : desc_) total_ += c;  // sum fits int64 -- up to ~3e14.
        // largest first, so desc_[q-1] is the q-th largest -- the bar q frees.
        std::sort(desc_.begin(), desc_.end(), std::greater<std::int64_t>());
    }

    // q in [2, n]. O(1) -- the sort already did the work.
    std::int64_t answer(int q) const {
        return total_ - desc_[static_cast<std::size_t>(q) - 1];
    }

    std::int64_t total() const { return total_; }

private:
    std::vector<std::int64_t> desc_;  // costs, largest first.
    std::int64_t total_;
};

// answer a batch of coupons in input order. sort once, O(1) per coupon.
inline std::vector<std::int64_t> solve(const std::vector<std::int64_t>& costs,
                                       const std::vector<int>& queries) {
    Discount d(costs);
    std::vector<std::int64_t> out;
    out.reserve(queries.size());
    for (int q : queries) out.push_back(d.answer(q));
    return out;
}

}  // namespace p0002
