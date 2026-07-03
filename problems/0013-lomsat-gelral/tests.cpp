#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

struct Tree {
    int n = 0;
    std::vector<int> color;                     // 0-indexed by vertex-1.
    std::vector<std::pair<int, int>> edges;     // 1-based undirected.
    bool operator==(const Tree&) const = default;
};

// map n random draws over [1, alpha] onto vertex colors -- small alphabets on
// purpose so counts tie and "sum every dominating color" bites.
std::vector<int> to_colors(int n, std::int64_t alpha, cpf::Rng& rng) {
    std::vector<int> c(static_cast<std::size_t>(n));
    for (auto& x : c) x = static_cast<int>(rng.in_range(1, alpha));
    return c;
}

// three shapes cover the degenerate ends of the sack recursion:
//   path    -- every vertex is its parent's heavy child, deepest possible chain.
//   star    -- root has n-1 light leaves, the heaviest single merge.
//   random  -- parent[i] uniform in [1, i-1], mixed light/heavy edges.
Tree make_path(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t;
    t.n = n;
    t.color = to_colors(n, alpha, rng);
    for (int v = 2; v <= n; ++v) t.edges.emplace_back(v - 1, v);
    return t;
}
Tree make_star(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t;
    t.n = n;
    t.color = to_colors(n, alpha, rng);
    for (int v = 2; v <= n; ++v) t.edges.emplace_back(1, v);
    return t;
}
Tree make_random(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t;
    t.n = n;
    t.color = to_colors(n, alpha, rng);
    for (int v = 2; v <= n; ++v)
        t.edges.emplace_back(static_cast<int>(rng.in_range(1, v - 1)), v);
    return t;
}

std::vector<std::int64_t> sol(const Tree& t) {
    return p0013::solve(t.n, t.color, t.edges);
}
std::vector<std::int64_t> ref(const Tree& t) {
    return p0013::reference_solve(t.n, t.color, t.edges);
}

}  // namespace

int main() {
    using namespace p0013;

    // ---- the statement example, verbatim. root 1, colors 1..4, a claw at 2. ----
    {
        std::vector<int> color = {1, 2, 3, 4};
        std::vector<std::pair<int, int>> edges = {{1, 2}, {2, 3}, {2, 4}};
        auto got = solve(4, color, edges);
        CPF_EQ(got.size(), std::size_t{4});
        CPF_EQ(got[0], std::int64_t{10});  // whole tree, all four colors once.
        CPF_EQ(got[1], std::int64_t{9});   // {2,3,4} once each.
        CPF_EQ(got[2], std::int64_t{3});   // leaf 3.
        CPF_EQ(got[3], std::int64_t{4});   // leaf 4.
        // the oracle must agree element for element -- ground truth, twice.
        CPF_CHECK(got == reference_solve(4, color, edges));
    }

    // ---- single vertex. its own color dominates, nothing else exists. ----
    {
        auto got = solve(1, {7}, {});
        CPF_EQ(got.size(), std::size_t{1});
        CPF_EQ(got[0], std::int64_t{7});
    }

    // ---- path of 5, all distinct colors. each subtree is a suffix chain; every
    //      color appears once so every subtree sums all its colors. ----
    {
        std::vector<int> color = {1, 2, 3, 4, 5};
        std::vector<std::pair<int, int>> edges = {{1, 2}, {2, 3}, {3, 4}, {4, 5}};
        auto got = solve(5, color, edges);
        // subtree(v) = {v, v+1, .., 5}, all distinct -> sum of those colors.
        CPF_EQ(got[0], std::int64_t{15});  // 1+2+3+4+5
        CPF_EQ(got[1], std::int64_t{14});  // 2+3+4+5
        CPF_EQ(got[2], std::int64_t{12});  // 3+4+5
        CPF_EQ(got[3], std::int64_t{9});   // 4+5
        CPF_EQ(got[4], std::int64_t{5});   // 5
        CPF_CHECK(got == reference_solve(5, color, edges));
    }

    // ---- star of 5, all same color. root sees color 1 five times, each leaf
    //      once -- a single color, so the sum is just that color everywhere. ----
    {
        std::vector<int> color = {1, 1, 1, 1, 1};
        std::vector<std::pair<int, int>> edges = {{1, 2}, {1, 3}, {1, 4}, {1, 5}};
        auto got = solve(5, color, edges);
        for (int i = 0; i < 5; ++i) CPF_EQ(got[static_cast<std::size_t>(i)], std::int64_t{1});
        CPF_CHECK(got == reference_solve(5, color, edges));
    }

    // ---- star with a forced tie at the root. two colors appear twice each among
    //      the leaves, the root once more of one -- exercises the tie sum. ----
    {
        // vertices 2..5 colored {2,2,3,3}; root colored 3. subtree(1): color 3
        // three times, color 2 twice -> only 3 dominates -> 3.
        std::vector<int> color = {3, 2, 2, 3, 3};
        std::vector<std::pair<int, int>> edges = {{1, 2}, {1, 3}, {1, 4}, {1, 5}};
        auto got = solve(5, color, edges);
        CPF_EQ(got[0], std::int64_t{3});
        CPF_CHECK(got == reference_solve(5, color, edges));
        // now flip root to a fresh color so 2 and 3 tie at two each -> 2+3 = 5.
        std::vector<int> color2 = {4, 2, 2, 3, 3};
        auto got2 = solve(5, color2, edges);
        CPF_EQ(got2[0], std::int64_t{5});
        CPF_CHECK(got2 == reference_solve(5, color2, edges));
    }

    // ---- int64 headroom: a fat subtree of distinct large colors. n distinct
    //      colors near n=1e5 sum past 2^31, so the answer must not be a 32-bit int.
    {
        const int n = 100000;
        std::vector<int> color(static_cast<std::size_t>(n));
        std::vector<std::pair<int, int>> edges;
        for (int v = 1; v <= n; ++v) color[static_cast<std::size_t>(v - 1)] = v;  // all distinct.
        for (int v = 2; v <= n; ++v) edges.emplace_back(1, v);  // star -> subtree(1) is everything.
        auto got = solve(n, color, edges);
        // subtree(1): every color once, so the answer is 1+2+..+n.
        std::int64_t want = static_cast<std::int64_t>(n) * (n + 1) / 2;  // 5000050000
        CPF_EQ(got[0], want);
        CPF_CHECK(want > (std::int64_t{1} << 31));  // proves 32-bit would wrap.
    }

    // ---- depth stress: a path at n=1e5. the solve is iterative, so a chain this
    //      deep must run to completion, not smash the call stack. ----
    {
        cpf::Rng rng(99);
        Tree t = make_path(100000, 4, rng);
        auto got = sol(t);
        CPF_EQ(got.size(), std::size_t{100000});
        // deepest leaf: subtree is one vertex, answer is its own color.
        CPF_EQ(got.back(), std::int64_t{t.color.back()});
    }

    // ---- randomized differential vs the O(n^2) oracle, thousands of trees. ----
    // small n and a tight alphabet so ties are the common case, not the exception;
    // one gen rotates through path / star / random shapes by seed.
    {
        const int kCases = 8000;
        bool ok = cpf::differential(
            kCases, 0x600Eu,
            [](cpf::Rng& rng) {
                const int n = static_cast<int>(rng.in_range(1, 60));
                const std::int64_t alpha = rng.in_range(1, 5);
                switch (rng.in_range(0, 2)) {
                    case 0: return make_path(n, alpha, rng);
                    case 1: return make_star(n, alpha, rng);
                    default: return make_random(n, alpha, rng);
                }
            },
            [](const Tree& t) { return sol(t); },
            [](const Tree& t) { return ref(t); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d trees, n<=60, alpha 1..5, "
                    "path/star/random\n",
                    kCases);
    }

    // ---- bigger random trees, wider alphabet -- fewer cases, more vertices, so
    //      the euler-block merge and heavy-path charge get a real workout. ----
    {
        const int kCases = 400;
        bool ok = cpf::differential(
            kCases, 0xBADCAFEu,
            [](cpf::Rng& rng) {
                const int n = static_cast<int>(rng.in_range(200, 500));
                const std::int64_t alpha = rng.in_range(2, 30);
                return make_random(n, alpha, rng);
            },
            [](const Tree& t) { return sol(t); },
            [](const Tree& t) { return ref(t); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d trees, n 200..500, alpha 2..30, "
                    "random\n",
                    kCases);
    }

    return cpf::report();
}
