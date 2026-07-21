#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0047::Edge;
using p0047::min_arborescence;
using p0047::min_arborescence_brute;

std::int64_t solve(int n, const std::vector<Edge>& e, int root) {
    return min_arborescence(n, e, root);
}

}  // namespace

int main() {
    // ---- hand examples ----

    // a triangle cycle that must contract. the two cheap edges 2<-3 and 3<-2
    // close a cycle, so greedy loops -- the only tree pays 10 to enter at 2.
    // cycle sum 1+1=2, entry reduced 10-1=9, total 11.
    CPF_EQ(solve(3, {{1, 2, 10}, {2, 3, 1}, {3, 2, 1}}, 1),
           static_cast<std::int64_t>(11));

    // a simple out-tree -- one legal shape, weights just add: 5+7+2.
    CPF_EQ(solve(4, {{1, 2, 5}, {2, 3, 7}, {1, 4, 2}}, 1),
           static_cast<std::int64_t>(14));

    // an unreachable vertex. 3 has no way in from the root -- no arborescence.
    CPF_EQ(solve(3, {{1, 2, 1}}, 1), static_cast<std::int64_t>(-1));

    // parallel edges -- the picker has to take the cheapest of the three.
    CPF_EQ(solve(2, {{1, 2, 9}, {1, 2, 3}, {1, 2, 7}}, 1),
           static_cast<std::int64_t>(3));

    // ---- edges ----

    // lone root, nothing to connect -- weight zero.
    CPF_EQ(solve(1, {}, 1), static_cast<std::int64_t>(0));

    // no edges but more than one vertex -- vertex 2 is stranded, -1.
    CPF_EQ(solve(2, {}, 1), static_cast<std::int64_t>(-1));

    // self-loops are noise. drop 2->2 and 1->1, the real answer is the 1->2 of 5.
    CPF_EQ(solve(2, {{2, 2, 1}, {1, 2, 5}, {2, 2, 100}, {1, 1, 3}}, 1),
           static_cast<std::int64_t>(5));

    // a non-1 root. root 2 reaches 3 and 1, cheapest tree 4+6.
    CPF_EQ(solve(3, {{2, 3, 4}, {3, 1, 6}, {1, 2, 100}}, 2),
           static_cast<std::int64_t>(10));

    // nested cycles: a whole 3-cycle 2->3->4->2 fed by one entry from the root.
    // cycle sum 1+1+1=3, entry 5 reduced by 1 -> 4, total 7.
    CPF_EQ(solve(4, {{1, 2, 5}, {2, 3, 1}, {3, 4, 1}, {4, 2, 1}}, 1),
           static_cast<std::int64_t>(7));

    // big weights near the cap -- the sum must not overflow an int.
    CPF_EQ(solve(3, {{1, 2, 1000000000LL}, {2, 3, 1000000000LL}}, 1),
           static_cast<std::int64_t>(2000000000LL));

    // ---- randomized differential vs the brute oracle ----
    //
    // thousands of small digraphs. weights are tiny on purpose -- ties and
    // 2-cycles show up constantly, which is exactly the contraction path. n stays
    // <= 7 so the in-edge-product oracle is affordable, and m is capped so the
    // oracle's enumeration can't blow up.
    {
        constexpr int kCases = 40000;
        int ran = 0, diffs = 0, neg_one = 0, contracted = 0;
        std::uint64_t first_bad = 0;

        for (int t = 0; t < kCases; ++t) {
            cpf::Rng rng(1000 + static_cast<std::uint64_t>(t));
            const int n = static_cast<int>(rng.in_range(1, 7));
            const int root = static_cast<int>(rng.in_range(1, n));
            const int cap = 3 * n < 15 ? 3 * n : 15;
            const int m = static_cast<int>(rng.in_range(0, cap));

            std::vector<Edge> e;
            e.reserve(static_cast<std::size_t>(m));
            for (int i = 0; i < m; ++i) {
                const int u = static_cast<int>(rng.in_range(1, n));
                const int v = static_cast<int>(rng.in_range(1, n));
                const std::int64_t w = rng.in_range(1, 4);  // tiny -> ties/cycles.
                e.push_back(Edge{u, v, w});
            }

            const std::int64_t got = min_arborescence(n, e, root);
            const std::int64_t want = min_arborescence_brute(n, e, root);
            ++ran;
            ++cpf::state().total;
            if (got != want) {
                ++cpf::state().failed;
                if (diffs == 0) first_bad = 1000 + static_cast<std::uint64_t>(t);
                ++diffs;
            }
            if (want == -1) ++neg_one;

            // did the graph actually force a contraction? a min-in-edge cycle
            // among the non-root vertices is what the plain greedy can't handle --
            // count those so the run proves it exercised the repair path.
            {
                std::vector<std::int64_t> in(static_cast<std::size_t>(n + 1),
                                             (1LL << 62));
                std::vector<int> pre(static_cast<std::size_t>(n + 1), 0);
                for (const Edge& x : e) {
                    if (x.u != x.v && x.w < in[static_cast<std::size_t>(x.v)]) {
                        in[static_cast<std::size_t>(x.v)] = x.w;
                        pre[static_cast<std::size_t>(x.v)] = x.u;
                    }
                }
                bool cyc = false;
                for (int s = 1; s <= n && !cyc; ++s) {
                    if (s == root ||
                        in[static_cast<std::size_t>(s)] == (1LL << 62))
                        continue;
                    int x = s, steps = 0;
                    while (x != root && steps <= n) {
                        if (in[static_cast<std::size_t>(x)] == (1LL << 62)) break;
                        x = pre[static_cast<std::size_t>(x)];
                        ++steps;
                    }
                    if (x != root && steps > n) cyc = true;
                }
                if (cyc && want != -1) ++contracted;
            }
        }

        CPF_CHECK(diffs == 0);
        if (diffs)
            std::printf("  first diff at seed %llu\n",
                        static_cast<unsigned long long>(first_bad));
        std::printf(
            "differential: %d digraphs, %d diffs, %d were -1, %d forced a "
            "contraction\n",
            ran, diffs, neg_one, contracted);
    }

    return cpf::report();
}
