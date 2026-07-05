#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

using p0033::Edge;
using p0033::i64;
using p0033::TreeMatch;
using p0033::TreeMatchRef;
using Edges = std::vector<Edge>;

namespace {

// build a tree by picking each new vertex's parent among the earlier ones. the
// PARENT RULE decides the shape:
//   "random"     -> parent uniform in [1, i-1]      (bushy random trees)
//   "path"       -> parent = i-1                     (a bamboo)
//   "star"       -> parent = 1                       (one hub)
//   "caterpillar"-> even i attach to a moving spine, odd i hang off as legs
// weights are drawn in [wlo, whi] so matchings of different sizes and the whole
// concavity structure both get stressed.
enum class Shape { Random, Path, Star, Caterpillar };

Edges make_tree(int n, Shape shape, cpf::Rng& r, i64 wlo, i64 whi) {
    Edges e;
    if (n <= 1) return e;
    e.reserve(static_cast<std::size_t>(n - 1));
    int spine = 1;  // caterpillar's current spine tip.
    for (int i = 2; i <= n; ++i) {
        int p;
        switch (shape) {
            case Shape::Path: p = i - 1; break;
            case Shape::Star: p = 1; break;
            case Shape::Caterpillar:
                if (i % 2 == 0) { p = spine; spine = i; }  // extend the spine.
                else { p = spine; }                         // hang a leg.
                break;
            case Shape::Random:
            default: p = static_cast<int>(r.in_range(1, i - 1)); break;
        }
        const i64 w = r.in_range(wlo, whi);
        e.push_back({i, p, w});
    }
    return e;
}

// solve every K in [0, n/2] and diff against the exact oracle. returns the count
// of (tree, K) pairs checked; asserts zero disagreement.
int sweep_all_k(int n, const Edges& edges) {
    TreeMatch fast(n, edges);
    TreeMatchRef oracle(n, edges);
    const std::vector<i64> want = oracle.all();  // want[k] = exact g(k) or -1.
    int cases = 0;
    for (int K = 0; K <= n / 2; ++K) {
        ++cases;
        CPF_EQ(fast.solve(K), want[static_cast<std::size_t>(K)]);
    }
    return cases;
}

}  // namespace

int main() {
    // --- hand examples from the statement. ---
    // path a-b-c-d, edge weights 1,5,1. K=1 grabs the middle 5. K=2 must take the
    // two outer 1s (only disjoint pair) = 2. K=3 impossible on 4 vertices -> -1.
    {
        Edges e{{1, 2, 1}, {2, 3, 5}, {3, 4, 1}};
        CPF_EQ(p0033::solve(0, 4, e), i64{0});
        CPF_EQ(p0033::solve(1, 4, e), i64{5});
        CPF_EQ(p0033::solve(2, 4, e), i64{2});
        CPF_EQ(p0033::solve(3, 4, e), i64{-1});
    }

    // star: hub 1 with legs of weight 3,7,4. any matching is a single edge, so
    // K=1 is the heaviest leg = 7, K>=2 impossible, K=0 is 0.
    {
        Edges e{{1, 2, 3}, {1, 3, 7}, {1, 4, 4}};
        CPF_EQ(p0033::solve(0, 4, e), i64{0});
        CPF_EQ(p0033::solve(1, 4, e), i64{7});
        CPF_EQ(p0033::solve(2, 4, e), i64{-1});
        CPF_EQ(p0033::solve(3, 4, e), i64{-1});
    }

    // --- edges. ---
    // single vertex, no edges: K=0 is 0, any K>=1 is -1.
    CPF_EQ(p0033::solve(0, 1, Edges{}), i64{0});
    CPF_EQ(p0033::solve(1, 1, Edges{}), i64{-1});
    // n=2, one edge of weight 9: K=0 -> 0, K=1 -> 9, K=2 -> -1.
    {
        Edges e{{1, 2, 9}};
        CPF_EQ(p0033::solve(0, 2, e), i64{0});
        CPF_EQ(p0033::solve(1, 2, e), i64{9});
        CPF_EQ(p0033::solve(2, 2, e), i64{-1});
    }
    // all-equal weights on a path of 6: g(k) = 4*k (k edges each worth 4), max
    // matching 3. K=3 -> 12, K=4 -> -1.
    {
        Edges e{{1, 2, 4}, {2, 3, 4}, {3, 4, 4}, {4, 5, 4}, {5, 6, 4}};
        CPF_EQ(p0033::solve(1, 6, e), i64{4});
        CPF_EQ(p0033::solve(2, 6, e), i64{8});
        CPF_EQ(p0033::solve(3, 6, e), i64{12});
        CPF_EQ(p0033::solve(4, 6, e), i64{-1});
    }
    // caterpillar by hand: spine 1-2-3, legs 4 off 2 and 5 off 3.
    // edges: 1-2 (w2), 2-3 (w2), 2-4 (w10), 3-5 (w10). max matching 2.
    // K=1 -> heaviest edge 10. K=2 -> 2-4 and 3-5 = 20 (disjoint, both heavy).
    {
        Edges e{{1, 2, 2}, {2, 3, 2}, {2, 4, 10}, {3, 5, 10}};
        CPF_EQ(p0033::solve(1, 5, e), i64{10});
        CPF_EQ(p0033::solve(2, 5, e), i64{20});
        CPF_EQ(p0033::solve(3, 5, e), i64{-1});
    }

    // --- int64 headroom. a path of many heavy edges: 2m+1 vertices, edges of
    // weight 1e9 on odd positions and 1 on even. exactly-m matching grabs all the
    // heavy ones -> m*1e9, past 2^31 and climbing toward 1e14. ---
    {
        const int m = 60;
        const int n = 2 * m + 1;
        Edges e;
        e.reserve(static_cast<std::size_t>(n - 1));
        for (int i = 1; i < n; ++i) {
            const i64 w = (i % 2 == 1) ? 1000000000LL : 1;  // heavy, light, heavy...
            e.push_back({i, i + 1, w});
        }
        // the m disjoint heavy edges (1-2, 3-4, ...) are a valid m-matching.
        CPF_EQ(p0033::solve(m, n, e), static_cast<i64>(m) * 1000000000LL);
        TreeMatch fast(n, e);
        TreeMatchRef oracle(n, e);
        const std::vector<i64> want = oracle.all();
        for (int K = 0; K <= n / 2; ++K)
            CPF_EQ(fast.solve(K), want[static_cast<std::size_t>(K)]);
        std::printf("int64 heavy-path: swept K=0..%d\n", n / 2);
    }

    // --- explicit shape sweeps over every K, hand-built. ---
    {
        cpf::Rng r(0xC0FFEE01);
        int cases = 0;
        const Shape shapes[4] = {Shape::Random, Shape::Path, Shape::Star,
                                 Shape::Caterpillar};
        for (Shape s : shapes)
            for (int n = 1; n <= 20; ++n)
                cases += sweep_all_k(n, make_tree(n, s, r, 1, 9));
        std::printf("shape sweeps: %d (tree,K) cases\n", cases);
    }

    // --- the main event: aliens vs the exact O(n^2) tree-knapsack oracle, ALL K
    // from 0..n/2, on thousands of small trees across four shapes and four weight
    // regimes. tight weight ranges force plateaus and count jumps; the -1 boundary
    // gets hit whenever K clears the tree's max matching. every tree is swept over
    // every K. any disagreement is a bug in the ALIENS code, never the oracle. ---
    {
        int trees = 0, cases = 0;
        const Shape shapes[4] = {Shape::Random, Shape::Path, Shape::Star,
                                 Shape::Caterpillar};
        // weight regimes: wide (distinct slopes), tight (many ties/plateaus),
        // binary (lots of equal marginals), single (all identical).
        const i64 wr[4][2] = {{1, 1000000000LL}, {1, 3}, {1, 2}, {5, 5}};
        for (int shape = 0; shape < 4; ++shape) {
            for (int reg = 0; reg < 4; ++reg) {
                for (std::uint64_t seed = 0; seed < 500; ++seed) {
                    cpf::Rng r(0xA11E5C0DEULL + shape * 1000000ULL +
                               reg * 100000ULL + seed);
                    const int n = static_cast<int>(r.in_range(1, 16));
                    const Edges e = make_tree(n, shapes[shape], r, wr[reg][0],
                                              wr[reg][1]);
                    ++trees;
                    cases += sweep_all_k(n, e);
                }
            }
        }
        std::printf("aliens vs oracle, all K: %d trees, %d (tree,K) cases\n",
                    trees, cases);
    }

    // --- a few larger random trees, still oracle-checkable, to stress deeper
    // recursion-free folds and bigger matchings. ---
    {
        cpf::Rng r(0xD00D5EED);
        int trees = 0, cases = 0;
        for (std::uint64_t seed = 0; seed < 200; ++seed) {
            const int n = static_cast<int>(r.in_range(30, 60));
            const Shape s = static_cast<Shape>(r.in_range(0, 3));
            const Edges e = make_tree(n, s, r, 1, 1000);
            ++trees;
            cases += sweep_all_k(n, e);
        }
        std::printf("larger trees, all K: %d trees, %d (tree,K) cases\n", trees,
                    cases);
    }

    return cpf::report();
}
