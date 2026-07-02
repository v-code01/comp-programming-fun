#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using namespace p0002;

namespace {

// one random instance: a cost array plus the coupon sizes to probe.
struct Input {
    std::vector<std::int64_t> costs;
    std::vector<int> queries;
};

// run every coupon through the O(n^2) oracle, in query order.
std::vector<std::int64_t> ref_solve(const Input& in) {
    ReferenceSolver ref;
    std::vector<std::int64_t> out;
    out.reserve(in.queries.size());
    for (int q : in.queries) out.push_back(ref.solve(in.costs, q));
    return out;
}

// probe every rank -- coupons 2..n. distinctness is a statement rule, not a
// correctness one, since each coupon answers on its own; testing all of them is
// the strongest sweep of the order statistics.
std::vector<int> all_coupons(int n) {
    std::vector<int> qs;
    qs.reserve(static_cast<std::size_t>(n - 1));
    for (int q = 2; q <= n; ++q) qs.push_back(q);
    return qs;
}

}  // namespace

int main() {
    // statement example: n=7, a=7 1 3 1 4 10 8. total 34, sorted 10 8 7 4 3 1 1.
    // coupon 3 frees the 3rd largest (7) -> 27; coupon 4 frees the 4th (4) -> 30.
    {
        std::vector<std::int64_t> a{7, 1, 3, 1, 4, 10, 8};
        Discount d(a);
        CPF_EQ(d.total(), 34LL);
        CPF_EQ(d.answer(3), 27LL);
        CPF_EQ(d.answer(4), 30LL);
    }

    // n=2, the floor. coupon 2 frees the cheaper of the two.
    {
        Discount d(std::vector<std::int64_t>{5, 8});
        CPF_EQ(d.answer(2), 8LL);  // 13 total, frees the 5.
    }

    // q=2 frees the second largest, q=n frees the smallest -- the two extremes.
    {
        std::vector<std::int64_t> a{4, 2, 9, 1, 6};  // total 22, desc 9 6 4 2 1.
        Discount d(a);
        CPF_EQ(d.answer(2), 22LL - 6);  // 2nd largest.
        CPF_EQ(d.answer(5), 22LL - 1);  // n-th largest = the min.
        CPF_EQ(d.answer(3), 22LL - 4);
    }

    // all-equal costs: every coupon frees the same value, whatever q is.
    {
        std::vector<std::int64_t> a(6, 7);  // total 42.
        Discount d(a);
        CPF_EQ(d.answer(2), 35LL);
        CPF_EQ(d.answer(6), 35LL);
    }

    // min and max costs together -- 1 and 1e9. no overflow: total is n*1e9 max,
    // ~3e14, well inside int64.
    {
        std::vector<std::int64_t> a{1, 1000000000LL, 1, 1000000000LL};
        Discount d(a);
        CPF_EQ(d.total(), 2000000002LL);
        CPF_EQ(d.answer(2), 2000000002LL - 1000000000LL);  // 2nd largest is 1e9.
        CPF_EQ(d.answer(4), 2000000002LL - 1);             // min is 1.
    }

    // full-size all-max: n=3e5 of 1e9. total 3e14, answer(q) = total - 1e9. this
    // is where a 32-bit accumulator would silently wrap -- int64 holds it.
    {
        std::vector<std::int64_t> a(300000, 1000000000LL);
        Discount d(a);
        CPF_EQ(d.total(), 300000LL * 1000000000LL);
        CPF_EQ(d.answer(2), 300000LL * 1000000000LL - 1000000000LL);
    }

    // randomized differential vs the oracle. small values first -- range 1..20
    // over n up to 40 packs the array with ties, the case the order-statistic
    // argument has to survive. every coupon 2..n on every instance.
    bool ties = cpf::differential(
        4000, 0xC0FFEEULL,
        [](cpf::Rng& rng) {
            int n = static_cast<int>(rng.in_range(2, 40));
            Input in;
            in.costs = rng.vec(static_cast<std::size_t>(n), 1, 20);
            in.queries = all_coupons(n);
            return in;
        },
        [](const Input& in) { return solve(in.costs, in.queries); }, ref_solve);
    CPF_CHECK(ties);

    // then the full value range -- 1..1e9, near-unique costs, to exercise the
    // int64 sums and the large-magnitude differences.
    bool wide = cpf::differential(
        4000, 0x5EEDULL,
        [](cpf::Rng& rng) {
            int n = static_cast<int>(rng.in_range(2, 40));
            Input in;
            in.costs = rng.vec(static_cast<std::size_t>(n), 1, 1000000000);
            in.queries = all_coupons(n);
            return in;
        },
        [](const Input& in) { return solve(in.costs, in.queries); }, ref_solve);
    CPF_CHECK(wide);

    std::printf("differential: 8000 instances, every coupon 2..n, 0 diffs\n");

    return cpf::report();
}
