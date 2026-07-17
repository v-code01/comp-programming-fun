#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0045::Beats;
using p0045::Reference;

// one op in the statement's alphabet, indices already 0-based.
//   type 1: a[i] = min(a[i], x) on [l, r]
//   type 2: a[i] += x on [l, r]        (x may be negative)
//   type 3: max of a[l..r]             (x unused)
//   type 4: sum of a[l..r]             (x unused)
struct Op {
    int type;
    int l;
    int r;
    std::int64_t x;
};

// replay a program on anything with the beats API, collecting the type-3 and
// type-4 answers in order. templated so the tree and the oracle walk the exact
// same path -- the only difference that can surface is a wrong answer.
template <typename DS>
std::vector<std::int64_t> run(DS& ds, const std::vector<std::int64_t>& init,
                              const std::vector<Op>& ops) {
    ds.build(init);
    std::vector<std::int64_t> out;
    for (const Op& o : ops) {
        if (o.type == 1)
            ds.chmin(o.l, o.r, o.x);
        else if (o.type == 2)
            ds.add(o.l, o.r, o.x);
        else if (o.type == 3)
            out.push_back(ds.query_max(o.l, o.r));
        else
            out.push_back(ds.query_sum(o.l, o.r));
    }
    return out;
}

// tree vs oracle on one fixed program -- both must emit the same answer stream.
void expect_agree(const std::vector<std::int64_t>& init,
                  const std::vector<Op>& ops) {
    Beats tree;
    Reference ref;
    CPF_CHECK(run(tree, init, ops) == run(ref, init, ops));
}

}  // namespace

int main() {
    // ---- hand example: a chmin that beats ---------------------------------
    // [5,4,3,2,1]; chmin[0,4] by 3 -> [3,3,3,2,1]; sum 12, max 3. then add
    // [0,2] +2 -> [5,5,5,2,1]; chmin[0,4] by 4 beats the three 5's -> [4,4,4,2,1]
    // (only the max band moves, se=2 stays put); sum 15, max 4.
    {
        std::vector<std::int64_t> init{5, 4, 3, 2, 1};
        std::vector<Op> ops{
            {1, 0, 4, 3},  // chmin by 3
            {4, 0, 4, 0},  // sum -> 12
            {3, 0, 4, 0},  // max -> 3
            {2, 0, 2, 2},  // add +2 on [0,2]
            {1, 0, 4, 4},  // chmin by 4 (beat)
            {4, 0, 4, 0},  // sum -> 15
            {3, 0, 4, 0},  // max -> 4
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        std::vector<std::int64_t> want{12, 3, 15, 4};
        CPF_CHECK(got == want);
        expect_agree(init, ops);
    }

    // ---- hand example: add then chmin, interleaved queries ----------------
    // [1,2,3,4,5]; add[1,3] +10 -> [1,12,13,14,5]; max 14; chmin[0,4] by 12
    // -> [1,12,12,12,5]; sum 42; max 12.
    {
        std::vector<std::int64_t> init{1, 2, 3, 4, 5};
        std::vector<Op> ops{
            {2, 1, 3, 10},  // add +10 on [1,3]
            {3, 0, 4, 0},   // max -> 14
            {1, 0, 4, 12},  // chmin by 12
            {4, 0, 4, 0},   // sum -> 42
            {3, 0, 4, 0},   // max -> 12
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        std::vector<std::int64_t> want{14, 42, 12};
        CPF_CHECK(got == want);
        expect_agree(init, ops);
    }

    // ---- edge: chmin above the max is a no-op -----------------------------
    {
        std::vector<std::int64_t> init{5, 9, 7};
        std::vector<Op> ops{
            {1, 0, 2, 1000000000LL},  // x >= mx, the prune returns at the root
            {4, 0, 2, 0},             // sum untouched -> 21
            {3, 0, 2, 0},             // max untouched -> 9
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK((got == std::vector<std::int64_t>{21, 9}));
        expect_agree(init, ops);
    }

    // ---- edge: chmin below everything flattens the whole array ------------
    {
        std::vector<std::int64_t> init{5, 9, 7, 3, 8};
        std::vector<Op> ops{
            {1, 0, 4, 2},  // x below every value -> all become 2
            {4, 0, 4, 0},  // sum -> 10
            {3, 0, 4, 0},  // max -> 2
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK((got == std::vector<std::int64_t>{10, 2}));
        expect_agree(init, ops);
    }

    // ---- edge: single element -- every op path on n = 1 -------------------
    {
        std::vector<std::int64_t> init{7};
        std::vector<Op> ops{
            {3, 0, 0, 0},   // max -> 7
            {1, 0, 0, 3},   // chmin by 3 -> 3
            {2, 0, 0, 10},  // add +10 -> 13
            {4, 0, 0, 0},   // sum -> 13
            {1, 0, 0, 20},  // chmin by 20, no-op -> 13
            {3, 0, 0, 0},   // max -> 13
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK((got == std::vector<std::int64_t>{7, 13, 13}));
        expect_agree(init, ops);
    }

    // ---- edge: all-equal array (cmx == len at the root) -------------------
    // every element is the max, so se stays -inf and a chmin caps the entire
    // band at once. then a partial chmin splits it so se becomes real.
    {
        std::vector<std::int64_t> init{6, 6, 6, 6};
        std::vector<Op> ops{
            {1, 0, 3, 4},  // cap all four -> [4,4,4,4]
            {4, 0, 3, 0},  // sum -> 16
            {1, 0, 1, 2},  // partial chmin -> [2,2,4,4], se now real
            {4, 0, 3, 0},  // sum -> 12
            {3, 0, 3, 0},  // max -> 4
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK((got == std::vector<std::int64_t>{16, 12, 4}));
        expect_agree(init, ops);
    }

    // ---- edge: negative adds and negative values --------------------------
    // add can drive values below zero; se must ride the shift, not the sentinel.
    {
        std::vector<std::int64_t> init{-3, 0, 5, 2};
        std::vector<Op> ops{
            {2, 0, 3, -4},  // -> [-7,-4,1,-2]
            {3, 0, 3, 0},   // max -> 1
            {1, 0, 3, -5},  // chmin by -5 -> [-7,-5,-5,-5]
            {4, 0, 3, 0},   // sum -> -22
            {2, 1, 3, 3},   // -> [-7,-2,-2,-2]
            {3, 0, 3, 0},   // max -> -2
        };
        Beats tree;
        std::vector<std::int64_t> got = run(tree, init, ops);
        CPF_CHECK((got == std::vector<std::int64_t>{1, -22, -2}));
        expect_agree(init, ops);
    }

    // ---- randomized differential vs the oracle ----------------------------
    // thousands of independent programs. small n keeps the O(n)-per-op oracle
    // cheap; a NARROW value band makes ties on the max common, so chmin lands in
    // the se < x < mx beat case constantly and the se/cmx bookkeeping and the
    // add+chmin lazy composition get hammered. a slice of ops reach for the full
    // 1e9 range and negatives so the sign and magnitude paths run too.
    {
        constexpr int kSeqs = 12000;
        cpf::Rng rng(0x5EA75ULL);
        int diffs = 0;
        std::uint64_t first_bad = 0;
        long long total_ops = 0;

        for (int s = 0; s < kSeqs; ++s) {
            int n = static_cast<int>(rng.in_range(1, 40));
            // narrow band by default so max-ties (and thus beats) are dense.
            bool wide = rng.in_range(0, 3) == 0;  // 1 in 4: stretch the range.
            std::int64_t vlo = wide ? -1000000000LL : -6;
            std::int64_t vhi = wide ? 1000000000LL : 6;

            std::vector<std::int64_t> init(static_cast<std::size_t>(n));
            for (auto& v : init) v = rng.in_range(vlo, vhi);

            int m = static_cast<int>(rng.in_range(1, 60));
            std::vector<Op> ops;
            ops.reserve(static_cast<std::size_t>(m));
            for (int qi = 0; qi < m; ++qi) {
                int t = static_cast<int>(rng.in_range(1, 4));
                int i = static_cast<int>(rng.in_range(0, n - 1));
                int j = static_cast<int>(rng.in_range(0, n - 1));
                int lo = i < j ? i : j;
                int hi = i < j ? j : i;
                std::int64_t x = 0;
                if (t == 1) {  // chmin: mostly in-band so it bites, sometimes not
                    x = rng.in_range(0, 4) == 0 ? rng.in_range(vlo, vhi)
                                                : rng.in_range(-6, 8);
                } else if (t == 2) {  // add: small, both signs
                    x = rng.in_range(-6, 6);
                }
                ops.push_back({t, lo, hi, x});
            }
            total_ops += m;

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
        std::printf("differential: %d programs, %lld ops, %d diffs\n", kSeqs,
                    total_ops, diffs);
    }

    return cpf::report();
}
