#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <set>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0039::Op;

struct Input {
    int n;
    std::vector<Op> ops;
};

// run the fast solution on a program.
std::vector<int> sol(const Input& in) { return p0039::solve(in.n, in.ops); }

// run the brute oracle on the same program.
std::vector<int> ref(const Input& in) { return p0039::referenceSolve(in.n, in.ops); }

// solution and oracle must emit the identical answer stream.
void expect_agree(const Input& in) { CPF_CHECK(sol(in) == ref(in)); }

// ---- a validity-tracking generator -----------------------------------------
// the solver trusts its input: type-2 removes only ever hit a currently present
// edge, and there's never a second copy of a live edge. so the generator carries a
// present-edge mirror and reads it before every choice. adds and removes churn the
// same handful of edges over and over, so each edge's life shatters into many short
// intervals -- exactly the fragmentation the segtree decomposition has to survive.
Input gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 10));
    Input in;
    in.n = n;

    std::set<std::pair<int, int>> present;  // the live-edge mirror.
    auto norm = [](int u, int v) {
        return u < v ? std::make_pair(u, v) : std::make_pair(v, u);
    };

    int m = static_cast<int>(rng.in_range(1, 60));
    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));

        if (roll < 35) {
            // add -- sample distinct pairs until one isn't already present.
            for (int tries = 0; tries < 30; ++tries) {
                int a = static_cast<int>(rng.in_range(1, n));
                int b = static_cast<int>(rng.in_range(1, n));
                if (a == b) continue;
                auto e = norm(a, b);
                if (present.count(e)) continue;
                present.insert(e);
                in.ops.push_back({1, a, b});
                break;
            }
        } else if (roll < 60) {
            // remove -- only valid if something is live. drop a random present edge.
            if (!present.empty()) {
                int idx = static_cast<int>(
                    rng.in_range(0, static_cast<std::int64_t>(present.size()) - 1));
                auto it = present.begin();
                std::advance(it, idx);
                in.ops.push_back({2, it->first, it->second});
                present.erase(it);
            }
        } else {
            // query -- any pair, u may equal v.
            int a = static_cast<int>(rng.in_range(1, n));
            int b = static_cast<int>(rng.in_range(1, n));
            in.ops.push_back({3, a, b});
        }
    }
    return in;
}

}  // namespace

int main() {
    // ---- hand example: add, connect, remove, disconnect --------------------
    {
        Input in{2, {{1, 1, 2}, {3, 1, 2}, {2, 1, 2}, {3, 1, 2}}};
        CPF_CHECK((sol(in) == std::vector<int>{1, 0}));
        expect_agree(in);
    }

    // ---- hand example: a triangle, add then tear down ----------------------
    // 1-2, 2-3, 1-3. connected throughout, then remove edges and watch it split.
    {
        Input in{3,
                 {
                     {1, 1, 2},  // {1,2}
                     {1, 2, 3},  // {1,2,3}
                     {1, 1, 3},  // still one component, now with a cycle
                     {3, 1, 3},  // -> 1
                     {2, 1, 3},  // drop the cycle edge -- path 1-2-3 remains
                     {3, 1, 3},  // -> 1
                     {2, 2, 3},  // now {1,2} and {3}
                     {3, 1, 3},  // -> 0
                     {3, 1, 2},  // -> 1
                 }};
        CPF_CHECK((sol(in) == std::vector<int>{1, 1, 0, 1}));
        expect_agree(in);
    }

    // ---- hand example: re-add after remove ---------------------------------
    // the same edge lives, dies, lives again -- two disjoint intervals for one pair.
    {
        Input in{2,
                 {
                     {1, 1, 2},  // add
                     {2, 1, 2},  // remove
                     {3, 1, 2},  // -> 0
                     {1, 1, 2},  // re-add
                     {3, 1, 2},  // -> 1
                 }};
        CPF_CHECK((sol(in) == std::vector<int>{0, 1}));
        expect_agree(in);
    }

    // ---- edge: query before any edge exists --------------------------------
    {
        Input in{3, {{3, 1, 2}}};
        CPF_CHECK((sol(in) == std::vector<int>{0}));
        expect_agree(in);
    }

    // ---- edge: self-query, u == v, is always 1 -----------------------------
    {
        Input in{3, {{3, 2, 2}, {1, 1, 2}, {3, 1, 1}}};
        CPF_CHECK((sol(in) == std::vector<int>{1, 1}));
        expect_agree(in);
    }

    // ---- edge: a chain that connects, then breaks in the middle ------------
    // 1-2-3-4-5. ends connected, then cut 3-4 and the two halves fall apart.
    {
        Input in{5,
                 {
                     {1, 1, 2}, {1, 2, 3}, {1, 3, 4}, {1, 4, 5},
                     {3, 1, 5},  // -> 1, whole chain
                     {2, 3, 4},  // cut the middle
                     {3, 1, 5},  // -> 0
                     {3, 1, 3},  // -> 1, left half
                     {3, 4, 5},  // -> 1, right half
                     {3, 2, 4},  // -> 0, straddles the cut
                 }};
        CPF_CHECK((sol(in) == std::vector<int>{1, 0, 1, 1, 0}));
        expect_agree(in);
    }

    // ---- edge: an edge left open runs to the end of the timeline -----------
    // never removed -- the interval must close at q, not vanish.
    {
        Input in{2, {{1, 1, 2}, {3, 1, 2}}};
        CPF_CHECK((sol(in) == std::vector<int>{1}));
        expect_agree(in);
    }

    // ---- randomized differential vs the brute oracle -----------------------
    // thousands of valid programs. adds and removes hammer the same edges so
    // intervals fragment hard; queries land at every instant. one disagreement
    // prints the seed that replays it.
    {
        constexpr int kSeqs = 8000;
        bool ok = cpf::differential(
            kSeqs, 0x0039A11,
            [](cpf::Rng& rng) { return gen(rng); },
            [](const Input& in) { return sol(in); },
            [](const Input& in) { return ref(in); });
        CPF_CHECK(ok);
        std::printf("differential: %d valid op-programs vs brute oracle\n", kSeqs);
    }

    return cpf::report();
}
