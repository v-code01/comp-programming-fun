#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

int main() {
    using namespace p0050;
    ReferenceSolver ref;

    // hand examples from the statement. [5,4,3,2,1] traced by hand:
    //   push 5                         -> {5}          cost 0
    //   push 4, top 5>4  pay 1         -> {4,4}        cost 1
    //   push 3, top 4>3  pay 1         -> {4,3,3}      cost 2
    //   push 2, top 4>2  pay 2         -> {3,3,2,2}    cost 4
    //   push 1, top 3>1  pay 2         -> {3,2,2,1,1}  cost 6
    // flatten toward the median, total 6.
    CPF_EQ(solve({5, 4, 3, 2, 1}), static_cast<i64>(6));
    CPF_EQ(solve({1, 2, 3, 4, 5}), static_cast<i64>(0));  // already sorted.
    CPF_EQ(solve({1, 3, 2, 4}), static_cast<i64>(1));     // one dip to smooth.
    CPF_EQ(solve({7, 7, 7, 7}), static_cast<i64>(0));     // all equal.

    // edges.
    CPF_EQ(solve({}), static_cast<i64>(0));      // empty -- nothing to fix.
    CPF_EQ(solve({42}), static_cast<i64>(0));    // n=1 -- always non-decreasing.
    CPF_EQ(solve({2, 1}), static_cast<i64>(1));  // two elements, pull one down.
    CPF_EQ(solve({1, 2}), static_cast<i64>(0));  // two, already fine.

    // negatives. [-3,-5,-1] -> pull -3 down to -5: cost 2.
    CPF_EQ(solve({-3, -5, -1}), static_cast<i64>(2));
    CPF_EQ(solve({-3, -5, -1}), ref.solve({-3, -5, -1}));
    CPF_EQ(solve({-1, -2, -3}), ref.solve({-1, -2, -3}));

    // int64 sanity: a big reverse-sorted pair drops ~2e9, must not wrap int32.
    CPF_EQ(solve({1000000000, -1000000000}), static_cast<i64>(2000000000));

    // differential vs the O(n*V) DP oracle on thousands of small arrays. small
    // value ranges pack in ties -- and the tie path (top == a_i, no pay; top >
    // a_i, pay) is exactly where the heap replacement can go wrong.
    cpf::Rng rng(20260713);
    int cases = 0, diffs = 0;
    auto check = [&](const std::vector<int>& v) {
        ++cases;
        if (solve(v) != ref.solve(v)) {
            ++diffs;
            if (diffs <= 3) {
                std::printf("  DIFF got %lld want %lld  [",
                            static_cast<long long>(solve(v)),
                            static_cast<long long>(ref.solve(v)));
                for (int x : v) std::printf(" %d", x);
                std::printf(" ]\n");
            }
        }
    };

    // heavy ties: tiny values, tiny n -- the merge fires constantly.
    for (int t = 0; t < 5000; ++t) {
        int n = static_cast<int>(rng.in_range(1, 15));
        std::vector<int> v(static_cast<std::size_t>(n));
        for (auto& x : v) x = static_cast<int>(rng.in_range(-5, 5));
        check(v);
    }
    // wider values, still small enough for the O(n*V) oracle.
    for (int t = 0; t < 3000; ++t) {
        int n = static_cast<int>(rng.in_range(1, 25));
        std::vector<int> v(static_cast<std::size_t>(n));
        for (auto& x : v) x = static_cast<int>(rng.in_range(-40, 40));
        check(v);
    }
    // already sorted -- oracle must agree AND the answer must be exactly 0.
    for (int t = 0; t < 1000; ++t) {
        int n = static_cast<int>(rng.in_range(1, 30));
        std::vector<int> v(static_cast<std::size_t>(n));
        for (auto& x : v) x = static_cast<int>(rng.in_range(-20, 20));
        std::sort(v.begin(), v.end());
        check(v);
        CPF_EQ(solve(v), static_cast<i64>(0));
    }
    // reverse sorted -- the worst case for the merge, every step pays.
    for (int t = 0; t < 1000; ++t) {
        int n = static_cast<int>(rng.in_range(1, 30));
        std::vector<int> v(static_cast<std::size_t>(n));
        for (auto& x : v) x = static_cast<int>(rng.in_range(-20, 20));
        std::sort(v.rbegin(), v.rend());
        check(v);
    }
    // all-equal -- cost 0, no matter the value.
    for (int t = 0; t < 500; ++t) {
        int n = static_cast<int>(rng.in_range(1, 30));
        int val = static_cast<int>(rng.in_range(-100, 100));
        std::vector<int> v(static_cast<std::size_t>(n), val);
        check(v);
        CPF_EQ(solve(v), static_cast<i64>(0));
    }

    CPF_CHECK(diffs == 0);
    std::printf("differential vs O(n*V) DP oracle: %d cases, %d diffs\n", cases,
                diffs);

    return cpf::report();
}
