#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Edges = std::vector<std::pair<int, int>>;

struct Instance {
    int n;
    Edges edges;
};

// one random small digraph. edges drawn freely with replacement, so multi-edges,
// cycles, back-edges, and unreachable vertices all show up on their own -- and
// each is a place Lengauer-Tarjan could quietly be wrong. self-loops are drawn
// too; both sides drop them, so they test that the drop is consistent.
Instance random_digraph(cpf::Rng& rng, int max_n, int max_extra) {
    const int n = static_cast<int>(rng.in_range(1, max_n));
    Instance in;
    in.n = n;
    // edge count spread from sparse (often disconnected) to dense (real merges).
    const int m = static_cast<int>(rng.in_range(0, n + max_extra));
    for (int e = 0; e < m; ++e) {
        const int u = static_cast<int>(rng.in_range(1, n));
        const int v = static_cast<int>(rng.in_range(1, n));
        in.edges.emplace_back(u, v);
    }
    return in;
}

// a compact printer for a failing case -- the seed replays it, but seeing the
// shape saves a round trip.
void dump(const Instance& in, const std::vector<int>& got,
          const std::vector<int>& want) {
    std::printf("  n=%d edges:", in.n);
    for (const auto& e : in.edges) std::printf(" %d->%d", e.first, e.second);
    std::printf("\n  got :");
    for (int v = 1; v <= in.n; ++v) std::printf(" %d", got[static_cast<std::size_t>(v)]);
    std::printf("\n  want:");
    for (int v = 1; v <= in.n; ++v) std::printf(" %d", want[static_cast<std::size_t>(v)]);
    std::printf("\n");
}

}  // namespace

int main() {
    using namespace p0036;
    ReferenceSolver ref;

    // ---- hand shapes. every expected vector is cross-checked against the
    // removal oracle so no number is just asserted from my head. ----

    // chain 1->2->3->4: each vertex's only way in is through the one before it.
    // idoms 0,1,2,3.
    {
        Edges e = {{1, 2}, {2, 3}, {3, 4}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[1], 0);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_EQ(got[4], 3);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // diamond 1->2, 1->3, 2->4, 3->4. two ways into 4, so neither 2 nor 3
    // dominates it -- the root does. idom(4)=1, and idom(2)=idom(3)=1.
    {
        Edges e = {{1, 2}, {1, 3}, {2, 4}, {3, 4}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[1], 0);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 1);
        CPF_EQ(got[4], 1);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // a cycle in the middle: 1->2->3->2, 3->4. the 2<->3 loop doesn't create a
    // second way in -- everything past 1 still funnels 1->2->3->4.
    {
        Edges e = {{1, 2}, {2, 3}, {3, 2}, {3, 4}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[1], 0);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_EQ(got[4], 3);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // an unreachable vertex. 1->2->3, and 4 dangles with an edge 4->2 that the
    // root never uses. idom(4)=0 -- not reachable, no dominator.
    {
        Edges e = {{1, 2}, {2, 3}, {4, 2}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[1], 0);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_EQ(got[4], 0);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // a merge where one branch is longer: 1->2->4, 1->3, 3->4. still two ways
    // into 4, so idom(4)=1.
    {
        Edges e = {{1, 2}, {2, 4}, {1, 3}, {3, 4}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[4], 1);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // ---- edges the machinery has to survive. ----

    // single vertex, no edges. only the root, idom 0.
    {
        Edges e = {};
        std::vector<int> got = dominator_tree(1, e);
        CPF_EQ(got[1], 0);
        CPF_CHECK(got == ref.solve(1, e));
    }
    // several vertices, no edges at all. only the root is reachable.
    {
        Edges e = {};
        std::vector<int> got = dominator_tree(5, e);
        for (int v = 1; v <= 5; ++v) CPF_EQ(got[v], 0);
        CPF_CHECK(got == ref.solve(5, e));
    }
    // a self-loop on the root and on a child -- must be ignored, not counted as
    // a second way in.
    {
        Edges e = {{1, 1}, {1, 2}, {2, 2}, {2, 3}};
        std::vector<int> got = dominator_tree(3, e);
        CPF_EQ(got[1], 0);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_CHECK(got == ref.solve(3, e));
    }
    // multi-edges: 1->2 twice. still one way in, idom(2)=1.
    {
        Edges e = {{1, 2}, {1, 2}, {2, 3}, {2, 3}};
        std::vector<int> got = dominator_tree(3, e);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_CHECK(got == ref.solve(3, e));
    }
    // two disjoint parts. root's component is 1->2->3; 4->5 is its own island,
    // both unreachable. idom(4)=idom(5)=0.
    {
        Edges e = {{1, 2}, {2, 3}, {4, 5}};
        std::vector<int> got = dominator_tree(5, e);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_EQ(got[4], 0);
        CPF_EQ(got[5], 0);
        CPF_CHECK(got == ref.solve(5, e));
    }
    // root with a back-edge into it: 1->2->3->1. the cycle through the root
    // changes nothing -- idoms stay the chain.
    {
        Edges e = {{1, 2}, {2, 3}, {3, 1}};
        std::vector<int> got = dominator_tree(3, e);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 2);
        CPF_CHECK(got == ref.solve(3, e));
    }
    // the nasty LT shape: a vertex reached both directly from the root and
    // through a longer detour, so its semidominator and idom differ and the
    // two-pass fixup actually has to fire.
    // 1->2, 1->3, 3->4, 4->2. 2 is reachable straight from 1 AND via 3->4->2,
    // so idom(2)=1 even though a tree path runs through 3 and 4.
    {
        Edges e = {{1, 2}, {1, 3}, {3, 4}, {4, 2}};
        std::vector<int> got = dominator_tree(4, e);
        CPF_EQ(got[2], 1);
        CPF_EQ(got[3], 1);
        CPF_EQ(got[4], 3);
        CPF_CHECK(got == ref.solve(4, e));
    }

    // ---- randomized differential vs the removal oracle. thousands of small
    // digraphs; cycles, unreachable vertices, and merges all arrive for free. ----
    {
        const int kCases = 40000;
        int diffs = 0;
        bool first = true;
        cpf::Rng rng(0x0036u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_digraph(rng, 14, 6);
            std::vector<int> got = dominator_tree(in.n, in.edges);
            std::vector<int> want = ref.solve(in.n, in.edges);
            if (got != want) {
                ++diffs;
                if (first) {
                    std::printf("  DIFF at case %d\n", c);
                    dump(in, got, want);
                    first = false;
                }
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf("randomized differential vs removal oracle: %d cases, "
                    "n<=14, %d diffs\n",
                    kCases, diffs);
    }

    // ---- a denser sweep: tiny n, many extra edges, so real merges and thick
    // cycles dominate and the semidominator relay gets hammered. count how many
    // cases carry an unreachable vertex, to prove that path is exercised. ----
    {
        const int kCases = 40000;
        int diffs = 0;
        int with_unreachable = 0;
        cpf::Rng rng(0xD0360u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_digraph(rng, 9, 18);
            std::vector<int> got = dominator_tree(in.n, in.edges);
            std::vector<int> want = ref.solve(in.n, in.edges);
            if (got != want) ++diffs;
            bool any_unreached = false;
            for (int v = 2; v <= in.n; ++v)
                if (want[static_cast<std::size_t>(v)] == 0) any_unreached = true;
            if (any_unreached) ++with_unreachable;
        }
        CPF_CHECK(diffs == 0);
        std::printf("dense/tiny differential: %d cases, n<=9, %d had an "
                    "unreachable vertex, %d diffs\n",
                    kCases, with_unreachable, diffs);
    }

    return cpf::report();
}
