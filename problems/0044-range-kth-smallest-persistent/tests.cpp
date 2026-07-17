#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0044::PersistentKth;
using p0044::ReferenceKth;

// build a fresh tree over a, run one query. small helper -- the hand examples
// each want their own array, so rebuilding per case is fine and clear.
std::int64_t kth(const std::vector<std::int64_t>& a, int l, int r, int k) {
    PersistentKth t;
    t.build(a);
    return t.query(l, r, k);
}

}  // namespace

int main() {
    // ---- the statement example, sorted by hand ------------------------------
    // a = [1,5,2,6,3,7,4]. query (2,5,3): a[2..5] = {5,2,6,3}, sorted {2,3,5,6},
    // the 3rd is 5.
    {
        std::vector<std::int64_t> a{1, 5, 2, 6, 3, 7, 4};
        CPF_EQ(kth(a, 2, 5, 3), 5);
        CPF_EQ(kth(a, 2, 5, 1), 2);  // k=1 -> min of the range.
        CPF_EQ(kth(a, 2, 5, 4), 6);  // k=len -> max of the range.
        CPF_EQ(kth(a, 1, 7, 1), 1);  // whole array, min.
        CPF_EQ(kth(a, 1, 7, 7), 7);  // whole array, max.
        CPF_EQ(kth(a, 1, 7, 4), 4);  // whole array, median -> 4.
        CPF_EQ(kth(a, 4, 4, 1), 6);  // single element -> itself.
    }

    // ---- edges enumerated ---------------------------------------------------

    // n=1: the only legal query is (1,1,1) -> the lone element.
    {
        std::vector<std::int64_t> a{42};
        CPF_EQ(kth(a, 1, 1, 1), 42);
    }

    // all equal: every value collapses to one compressed leaf. any k, any range
    // returns that value -- the duplicate-heavy path through compression.
    {
        std::vector<std::int64_t> a{7, 7, 7, 7, 7};
        CPF_EQ(kth(a, 1, 5, 1), 7);
        CPF_EQ(kth(a, 1, 5, 5), 7);
        CPF_EQ(kth(a, 2, 4, 2), 7);
    }

    // all distinct, ascending: k-th smallest of a[l..r] is just a[l+k-1].
    {
        std::vector<std::int64_t> a{10, 20, 30, 40, 50};
        CPF_EQ(kth(a, 1, 5, 3), 30);
        CPF_EQ(kth(a, 2, 4, 1), 20);
        CPF_EQ(kth(a, 2, 4, 3), 40);
    }

    // negatives, and a mix straddling zero -- signed values must compress and
    // order correctly. a = [-5,3,-5,0,-1], a[1..5] sorted = {-5,-5,-1,0,3}.
    {
        std::vector<std::int64_t> a{-5, 3, -5, 0, -1};
        CPF_EQ(kth(a, 1, 5, 1), -5);
        CPF_EQ(kth(a, 1, 5, 2), -5);  // the duplicate -5.
        CPF_EQ(kth(a, 1, 5, 3), -1);
        CPF_EQ(kth(a, 1, 5, 4), 0);
        CPF_EQ(kth(a, 1, 5, 5), 3);
        CPF_EQ(kth(a, 2, 4, 2), 0);   // {3,-5,0} sorted {-5,0,3}, 2nd -> 0.
    }

    // full-magnitude values -- 1e9 and -1e9 sit at the ends, no overflow, no
    // int32 surprise on the value side.
    {
        std::vector<std::int64_t> a{1000000000, -1000000000, 0};
        CPF_EQ(kth(a, 1, 3, 1), -1000000000);
        CPF_EQ(kth(a, 1, 3, 2), 0);
        CPF_EQ(kth(a, 1, 3, 3), 1000000000);
    }

    // ---- randomized differential vs the sort oracle -------------------------
    // thousands of small random arrays, each hit with a batch of random valid
    // queries. the value pool is deliberately narrow so duplicates are common,
    // and it dips negative so the sign path gets worked. every query hits k=1,
    // k=len, and single-element ranges by construction of the random spread. any
    // disagreement prints the seed that replays it.
    {
        constexpr int kArrays = 5000;
        cpf::Rng rng(0x0044ULL);
        int queries = 0;
        int diffs = 0;
        std::uint64_t first_bad = 0;

        for (int s = 0; s < kArrays; ++s) {
            const int n = static_cast<int>(rng.in_range(1, 60));
            // narrow, signed pool -> lots of ties and negatives.
            const std::int64_t lo = rng.in_range(-8, 0);
            const std::int64_t hi = lo + rng.in_range(0, 12);
            std::vector<std::int64_t> a(static_cast<std::size_t>(n));
            for (auto& x : a) x = rng.in_range(lo, hi);

            PersistentKth fast;
            fast.build(a);
            ReferenceKth ref;
            ref.build(a);

            const int m = static_cast<int>(rng.in_range(1, 40));
            for (int qi = 0; qi < m; ++qi) {
                const int l = static_cast<int>(rng.in_range(1, n));
                const int r = static_cast<int>(rng.in_range(l, n));
                const int k = static_cast<int>(rng.in_range(1, r - l + 1));
                ++queries;
                if (fast.query(l, r, k) != ref.query(l, r, k)) {
                    if (diffs == 0) first_bad = static_cast<std::uint64_t>(s);
                    ++diffs;
                }
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first differing array index: %llu\n",
                        static_cast<unsigned long long>(first_bad));
        }
        std::printf("differential: %d arrays, %d queries, %d diffs\n", kArrays,
                    queries, diffs);
    }

    // ---- larger arrays, wider domain, still vs the oracle --------------------
    // push n up so the tree grows several levels deep and the compression maps a
    // wide value range; the oracle stays O(len log len) per query so keep the
    // query batch modest.
    {
        constexpr int kArrays = 300;
        cpf::Rng rng(0xC0FFEEULL);
        int queries = 0;
        int diffs = 0;

        for (int s = 0; s < kArrays; ++s) {
            const int n = static_cast<int>(rng.in_range(200, 1200));
            std::vector<std::int64_t> a(static_cast<std::size_t>(n));
            for (auto& x : a) x = rng.in_range(-1000000000, 1000000000);

            PersistentKth fast;
            fast.build(a);
            ReferenceKth ref;
            ref.build(a);

            for (int qi = 0; qi < 50; ++qi) {
                const int l = static_cast<int>(rng.in_range(1, n));
                const int r = static_cast<int>(rng.in_range(l, n));
                const int k = static_cast<int>(rng.in_range(1, r - l + 1));
                ++queries;
                if (fast.query(l, r, k) != ref.query(l, r, k)) ++diffs;
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf("larger-array differential: %d arrays, %d queries, %d diffs\n",
                    kArrays, queries, diffs);
    }

    return cpf::report();
}
