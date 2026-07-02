#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0009::Beats;
using p0009::Reference;

// one op in the statement's alphabet, indices already 0-based.
//   type 1: sum a[l..r]      (x unused)
//   type 2: a[i] %= x on [l, r]
//   type 3: a[l] = x         (r unused)
struct Op {
    int type;
    int l;
    int r;
    std::int64_t x;
};

// replay a program on any structure with the beats API, collecting the type-1
// answers in order. templated so the tree and the oracle run the identical path.
template <typename DS>
std::vector<std::int64_t> run(DS& ds, const std::vector<std::int64_t>& init,
                              const std::vector<Op>& ops) {
    ds.build(init);
    std::vector<std::int64_t> out;
    for (const Op& o : ops) {
        if (o.type == 1)
            out.push_back(ds.range_sum(o.l, o.r));
        else if (o.type == 2)
            ds.range_mod(o.l, o.r, o.x);
        else
            ds.point_set(o.l, o.x);
    }
    return out;
}

// tree vs oracle on one fixed program -- both must emit the same type-1 stream.
void expect_agree(const std::vector<std::int64_t>& init,
                  const std::vector<Op>& ops) {
    Beats tree;
    Reference ref;
    CPF_CHECK(run(tree, init, ops) == run(ref, init, ops));
}

}  // namespace

int main() {
    // ---- statement example 1 (indices shifted 1 -> 0) ---------------------
    // a = 1 2 3 4 5 ; mod [3,5] by 4 ; set a3=5 ; sum[2,5]=8 ; mod[1,3] by 3 ;
    // sum[1,3]=5.
    {
        std::vector<std::int64_t> init{1, 2, 3, 4, 5};
        std::vector<Op> ops{
            {2, 2, 4, 4},  // mod [3,5] by 4
            {3, 2, 0, 5},  // set a[3] = 5
            {1, 1, 4, 0},  // sum [2,5] -> 8
            {2, 0, 2, 3},  // mod [1,3] by 3
            {1, 0, 2, 0},  // sum [1,3] -> 5
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        std::vector<std::int64_t> want{8, 5};
        CPF_CHECK(got == want);
        expect_agree(init, ops);
    }

    // ---- statement example 2 ----------------------------------------------
    {
        std::vector<std::int64_t> init{6, 9, 6, 7, 6, 1, 10, 10, 9, 5};
        std::vector<Op> ops{
            {1, 2, 8, 0},   // sum [3,9]  -> 49
            {2, 6, 9, 9},   // mod [7,10] by 9
            {2, 4, 9, 8},   // mod [5,10] by 8
            {1, 3, 6, 0},   // sum [4,7]  -> 15
            {3, 2, 0, 7},   // set a[3] = 7
            {2, 6, 8, 9},   // mod [7,9] by 9
            {1, 1, 3, 0},   // sum [2,4]  -> 23
            {1, 5, 5, 0},   // sum [6,6]  -> 1
            {1, 4, 8, 0},   // sum [5,9]  -> 9
            {3, 0, 0, 10},  // set a[1] = 10
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        std::vector<std::int64_t> want{49, 15, 23, 1, 9};
        CPF_CHECK(got == want);
        expect_agree(init, ops);
    }

    // ---- enumerated edges --------------------------------------------------
    // single element: mod reduces it, then a set repopulates it to the ceiling.
    {
        std::vector<std::int64_t> init{7};
        std::vector<Op> ops{
            {1, 0, 0, 0},              // sum -> 7
            {2, 0, 0, 3, },            // 7 %= 3 -> 1
            {1, 0, 0, 0},              // sum -> 1
            {3, 0, 0, 1000000000LL},   // set -> 1e9
            {1, 0, 0, 0},              // sum -> 1e9
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        std::vector<std::int64_t> want{7, 1, 1000000000LL};
        CPF_CHECK(got == want);
        expect_agree(init, ops);
    }

    // mod by 1: everything collapses to 0.
    {
        std::vector<std::int64_t> init{5, 9, 1000000000LL};
        expect_agree(init, {{2, 0, 2, 1}, {1, 0, 2, 0}});  // -> 0
        Beats tree;
        CPF_CHECK((run(tree, init, {{2, 0, 2, 1}, {1, 0, 2, 0}}) ==
                   std::vector<std::int64_t>{0}));
    }

    // mod by a huge x, larger than every value: the prune fires, no-op.
    {
        std::vector<std::int64_t> init{5, 9, 7};
        Beats tree;
        std::vector<std::int64_t> got =
            run(tree, init, {{2, 0, 2, 1000000000LL}, {1, 0, 2, 0}});
        CPF_CHECK(got == std::vector<std::int64_t>{21});  // untouched sum.
    }

    // point-set then mod: the set lifts a value back above x so the mod bites.
    {
        std::vector<std::int64_t> init{1, 1, 1};
        std::vector<Op> ops{
            {3, 1, 0, 100},  // a[1] = 100
            {2, 0, 2, 7},    // [1,100,1] %= 7 -> [1,2,1]
            {1, 0, 2, 0},    // sum -> 4
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK(got == std::vector<std::int64_t>{4});
        expect_agree(init, ops);
    }

    // all-equal array: mod by the value zeroes it; mod by a divisor shifts all.
    {
        std::vector<std::int64_t> init{5, 5, 5, 5};
        expect_agree(init, {{2, 0, 3, 5}, {1, 0, 3, 0}});  // -> 0
        expect_agree(init, {{2, 0, 3, 3}, {1, 0, 3, 0}});  // -> 8 (all become 2)
        Beats a;
        CPF_CHECK((run(a, init, {{2, 0, 3, 3}, {1, 0, 3, 0}}) ==
                   std::vector<std::int64_t>{8}));
    }

    // ---- randomized differential vs the oracle ----------------------------
    // thousands of independent programs. small n so the O(n)-per-op oracle stays
    // cheap; a value pool skewed small so mods actually fire, salted with the
    // full 1e9 range so the int64 sum path and the prune both get exercised.
    {
        constexpr int kSeqs = 8000;
        cpf::Rng rng(0xC1FDULL);
        int diffs = 0;
        std::uint64_t first_bad = 0;

        for (int s = 0; s < kSeqs; ++s) {
            int n = static_cast<int>(rng.in_range(1, 40));
            bool big_vals = rng.in_range(0, 3) == 0;  // 1 in 4: allow huge cells.
            std::int64_t vhi = big_vals ? 1000000000LL : 25;

            std::vector<std::int64_t> init(static_cast<std::size_t>(n));
            for (auto& v : init) v = rng.in_range(1, vhi);

            int m = static_cast<int>(rng.in_range(1, 60));
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int q = 0; q < m; ++q) {
                int t = static_cast<int>(rng.in_range(1, 3));
                int i = static_cast<int>(rng.in_range(0, n - 1));
                int j = static_cast<int>(rng.in_range(0, n - 1));
                int lo = i < j ? i : j;
                int hi = i < j ? j : i;
                if (t == 1) {
                    ops.push_back({1, lo, hi, 0});
                } else if (t == 2) {
                    // mostly small x so reductions land; sometimes big / =1.
                    std::int64_t x =
                        rng.in_range(0, 4) == 0 ? rng.in_range(1, vhi)
                                                : rng.in_range(1, 12);
                    ops.push_back({2, lo, hi, x});
                } else {
                    std::int64_t x = rng.in_range(1, vhi);
                    ops.push_back({3, lo, 0, x});  // set a[lo] = x
                }
            }

            Beats tree;
            Reference ref;
            if (run(tree, init, ops) != run(ref, init, ops)) {
                if (diffs == 0) first_bad = static_cast<std::uint64_t>(s);
                ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs)
            std::printf("  first differing sequence index: %llu\n",
                        static_cast<unsigned long long>(first_bad));
        std::printf("differential: %d programs, %d diffs\n", kSeqs, diffs);
    }

    return cpf::report();
}
