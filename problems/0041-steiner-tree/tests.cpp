#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0041::Edge;

struct Instance {
    int n = 0;
    std::vector<Edge> edges;
    std::vector<int> terminals;
};

// weight of the MST of the subgraph induced on the TERMINALS ALONE, or kInf if
// that subgraph is disconnected. this is the "no Steiner points allowed" answer.
// the Steiner tree is always <= it, and when it is STRICTLY less, a non-terminal
// vertex is doing real work -- which is the case the whole DP exists for. the
// sweeps below count those, so "we tested real Steiner points" is a number, not
// a claim.
std::int64_t terminals_only_mst(const Instance& in) {
    unsigned S = 0;
    for (int t : in.terminals) S |= 1u << t;
    std::vector<Edge> sorted = in.edges;
    std::sort(sorted.begin(), sorted.end(),
              [](const Edge& a, const Edge& b) { return a.w < b.w; });

    std::vector<int> parent(static_cast<std::size_t>(in.n));
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&parent](int x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(
                parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };

    int comps = static_cast<int>(in.terminals.size());
    std::int64_t total = 0;
    for (const Edge& e : sorted) {
        if (!((S >> e.u) & 1u) || !((S >> e.v) & 1u)) continue;
        const int ru = find(e.u);
        const int rv = find(e.v);
        if (ru == rv) continue;
        parent[static_cast<std::size_t>(ru)] = rv;
        total += e.w;
        --comps;
    }
    return comps == 1 ? total : p0041::kInf;
}

// a random CONNECTED weighted graph plus a random terminal set.
//
// spine first: vertex v attaches to a uniform earlier vertex, so the graph is
// connected by construction and every extra edge is optional. then sprinkle
// extras at the requested density. weights are small on purpose -- a wide weight
// range makes the shortest route obvious and the DP trivially right; tight
// weights make hub-vs-direct genuinely close, which is where the merge/grow
// ordering breaks if it's wrong.
Instance random_instance(cpf::Rng& rng, int max_n, int max_k, int extra_pct,
                         std::int64_t wlo, std::int64_t whi) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(1, max_n));

    for (int v = 1; v < in.n; ++v) {
        const int p = static_cast<int>(rng.in_range(0, v - 1));
        in.edges.push_back({p, v, rng.in_range(wlo, whi)});
    }
    for (int u = 0; u < in.n; ++u) {
        for (int v = u + 1; v < in.n; ++v) {
            if (rng.in_range(1, 100) <= extra_pct) {
                in.edges.push_back({u, v, rng.in_range(wlo, whi)});
            }
        }
    }

    std::vector<int> perm(static_cast<std::size_t>(in.n));
    std::iota(perm.begin(), perm.end(), 0);
    for (int i = in.n - 1; i > 0; --i) {  // fisher-yates on the seeded rng.
        const int j = static_cast<int>(rng.in_range(0, i));
        std::swap(perm[static_cast<std::size_t>(i)],
                  perm[static_cast<std::size_t>(j)]);
    }
    const int k = static_cast<int>(rng.in_range(1, std::min(max_k, in.n)));
    in.terminals.assign(perm.begin(), perm.begin() + k);
    return in;
}

// one sweep. runs the DP against the 2^n + MST oracle, and reports how many of
// the cases actually required a Steiner point.
void sweep(const char* label, int cases, std::uint64_t seed, int max_n,
           int max_k, int extra_pct, std::int64_t wlo, std::int64_t whi) {
    cpf::Rng rng(seed);
    int diffs = 0;
    int steiner_needed = 0;
    int first_bad = -1;
    for (int c = 0; c < cases; ++c) {
        const Instance in = random_instance(rng, max_n, max_k, extra_pct, wlo, whi);
        const std::int64_t got = p0041::steiner_tree(in.n, in.edges, in.terminals);
        const std::int64_t want =
            p0041::steiner_tree_brute(in.n, in.edges, in.terminals);
        if (got != want) {
            if (diffs == 0) first_bad = c;
            ++diffs;
        }
        if (in.terminals.size() >= 2 && got < terminals_only_mst(in)) {
            ++steiner_needed;
        }
    }
    ++cpf::state().total;
    if (diffs != 0) {
        ++cpf::state().failed;
        std::printf("  FAIL %s: %d diffs, first at case %d\n", label, diffs,
                    first_bad);
    }
    std::printf("%-34s %6d cases, %6d needed a Steiner point, %d diffs\n", label,
                cases, steiner_needed, diffs);
}

}  // namespace

int main() {
    using p0041::steiner_tree;
    using p0041::steiner_tree_brute;

    // ---- hand shapes. every expected number is also run through the oracle, so
    // nothing here is asserted purely from my head. ----

    // THE STAR. center 0 is not a terminal and there is no other route between
    // the leaves -- the center is a REQUIRED Steiner point. 1+1+1 = 3.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}};
        std::vector<int> t = {1, 2, 3};
        CPF_EQ(steiner_tree(n, e, t), 3);
        CPF_EQ(steiner_tree_brute(n, e, t), 3);
    }
    // same star, now the leaves ARE directly connected -- but expensively. the
    // terminal-only MST is 10 + 10 = 20; going through the center is 3. the
    // Steiner point wins on price, not on necessity.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 1},  {0, 2, 1},  {0, 3, 1},
                               {1, 2, 10}, {2, 3, 10}, {1, 3, 10}};
        std::vector<int> t = {1, 2, 3};
        CPF_EQ(steiner_tree(n, e, t), 3);
        CPF_EQ(steiner_tree_brute(n, e, t), 3);
    }
    // flip the prices: the hub is expensive, direct edges are cheap. now the
    // answer is the terminal-only MST (1 + 1 = 2) and the center is dead weight.
    // if the DP "always routes through the hub", this catches it.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 10}, {0, 2, 10}, {0, 3, 10},
                               {1, 2, 1},  {2, 3, 1},  {1, 3, 5}};
        std::vector<int> t = {1, 2, 3};
        CPF_EQ(steiner_tree(n, e, t), 2);
        CPF_EQ(steiner_tree_brute(n, e, t), 2);
    }

    // A PATH 0-1-2-3-4, terminals at both ends. everything between is a forced
    // Steiner point and the answer is the whole path: 4.
    {
        const int n = 5;
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
        std::vector<int> t = {0, 4};
        CPF_EQ(steiner_tree(n, e, t), 4);
        CPF_EQ(steiner_tree_brute(n, e, t), 4);
    }

    // k = 1 -> 0. a single terminal is already a connected subgraph.
    {
        const int n = 5;
        std::vector<Edge> e = {{0, 1, 7}, {1, 2, 7}, {2, 3, 7}, {3, 4, 7}};
        std::vector<int> t = {2};
        CPF_EQ(steiner_tree(n, e, t), 0);
        CPF_EQ(steiner_tree_brute(n, e, t), 0);
    }
    // n = 1, m = 0, k = 1 -- the smallest legal input. 0.
    {
        std::vector<Edge> e;
        std::vector<int> t = {0};
        CPF_EQ(steiner_tree(1, e, t), 0);
        CPF_EQ(steiner_tree_brute(1, e, t), 0);
    }
    // k = 0 is outside the stated constraints but the guard is one branch and a
    // caller with an empty terminal list deserves 0, not a read past the end.
    {
        std::vector<Edge> e = {{0, 1, 3}};
        std::vector<int> t;
        CPF_EQ(steiner_tree(2, e, t), 0);
        CPF_EQ(steiner_tree_brute(2, e, t), 0);
    }

    // TWO TERMINALS -> exactly the shortest path. the detour 0-1-2-3 costs 3, the
    // direct edge 0-3 costs 5. the DP with k = 2 is just Dijkstra wearing a hat.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {0, 3, 5}};
        std::vector<int> t = {0, 3};
        CPF_EQ(steiner_tree(n, e, t), 3);
        CPF_EQ(steiner_tree_brute(n, e, t), 3);
    }
    // two terminals, direct edge is cheaper than the long way around.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {0, 3, 5}};
        std::vector<int> t = {0, 3};
        CPF_EQ(steiner_tree(n, e, t), 5);
        CPF_EQ(steiner_tree_brute(n, e, t), 5);
    }

    // ALL VERTICES ARE TERMINALS -> the answer is the MST of the whole graph.
    // no vertex is free to drop, so the cheapest connected spanning subgraph is
    // the MST. this cycle's MST drops the single 9-edge: 1+2+3 = 6.
    {
        const int n = 4;
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 2}, {2, 3, 3}, {3, 0, 9}};
        std::vector<int> t = {0, 1, 2, 3};
        CPF_EQ(steiner_tree(n, e, t), 6);
        CPF_EQ(steiner_tree_brute(n, e, t), 6);
    }

    // A STEINER POINT THAT STRICTLY HELPS, with the hub cost split three ways.
    // four terminals 1..4 on a hub 0 at cost 3 each -> 12. the direct ring
    // between the terminals costs 5 an edge -> a terminal-only tree needs 3 of
    // them = 15. the DP has to find 12, which means finding the hub.
    {
        const int n = 5;
        std::vector<Edge> e = {{0, 1, 3}, {0, 2, 3}, {0, 3, 3}, {0, 4, 3},
                               {1, 2, 5}, {2, 3, 5}, {3, 4, 5}, {4, 1, 5}};
        std::vector<int> t = {1, 2, 3, 4};
        CPF_EQ(steiner_tree(n, e, t), 12);
        CPF_EQ(steiner_tree_brute(n, e, t), 12);
        // and the hub is doing real work -- the terminal-only tree really is 15.
        Instance in{n, e, t};
        CPF_EQ(terminals_only_mst(in), 15);
    }

    // TWO HUBS. terminals 2,3 sit under hub 0; terminals 4,5 under hub 1; the
    // hubs are joined by a bridge. the optimal tree uses BOTH Steiner points and
    // the merge transition has to meet at a hub, not at a terminal:
    // 1+1 (hub0's leaves) + 1 (bridge) + 1+1 (hub1's leaves) = 5.
    {
        const int n = 6;
        std::vector<Edge> e = {{0, 2, 1}, {0, 3, 1}, {0, 1, 1},
                               {1, 4, 1}, {1, 5, 1}};
        std::vector<int> t = {2, 3, 4, 5};
        CPF_EQ(steiner_tree(n, e, t), 5);
        CPF_EQ(steiner_tree_brute(n, e, t), 5);
    }

    // parallel edges and a self-loop -- both legal to hand a solver, both must be
    // ignored gracefully. the cheap parallel copy is the one that counts.
    {
        const int n = 3;
        std::vector<Edge> e = {{0, 1, 9}, {0, 1, 2}, {1, 2, 4}, {2, 2, 1}};
        std::vector<int> t = {0, 2};
        CPF_EQ(steiner_tree(n, e, t), 6);
        CPF_EQ(steiner_tree_brute(n, e, t), 6);
    }

    // WEIGHTS AT THE CEILING. a 100-vertex path at w = 1e6, terminals at the two
    // ends: 99 * 1e6 = 99'000'000. int32 survives this one, but the sum is the
    // thing to keep an eye on, and the DP's kInf arithmetic is what actually
    // needs the 64 bits.
    {
        const int n = 100;
        std::vector<Edge> e;
        for (int i = 0; i + 1 < n; ++i) e.push_back({i, i + 1, 1000000});
        std::vector<int> t = {0, n - 1};
        CPF_EQ(steiner_tree(n, e, t), 99000000LL);
    }
    // ten terminals spread along that same path -- the tree is the span between
    // the extreme two, and every vertex between them is a Steiner point.
    {
        const int n = 100;
        std::vector<Edge> e;
        for (int i = 0; i + 1 < n; ++i) e.push_back({i, i + 1, 1000000});
        std::vector<int> t = {5, 12, 19, 27, 33, 41, 58, 66, 74, 91};
        CPF_EQ(steiner_tree(n, e, t), (91 - 5) * 1000000LL);
    }

    // ---- randomized differential vs the 2^n + MST oracle. thousands of small
    // graphs across densities: sparse (trees and near-trees) is where forced
    // Steiner points live, dense is where the merge has real choices to make. ----
    sweep("sparse n<=10 k<=4 w<=10", 4000, 0x41A1u, 10, 4, 5, 1, 10);
    sweep("medium n<=10 k<=4 w<=10", 4000, 0x41A2u, 10, 4, 30, 1, 10);
    sweep("dense  n<=10 k<=4 w<=10", 4000, 0x41A3u, 10, 4, 70, 1, 10);
    sweep("tight  n<=12 k<=4 w<=3", 4000, 0x41A4u, 12, 4, 20, 1, 3);
    sweep("wide   n<=12 k<=4 w<=1e6", 2000, 0x41A5u, 12, 4, 25, 1, 1000000);
    sweep("tiny   n<=6  k<=4 w<=5", 4000, 0x41A6u, 6, 4, 40, 1, 5);
    // k up to 6 on n <= 12: pushes the merge tree deeper while the oracle stays
    // cheap. this is the sweep that actually stresses submask enumeration.
    sweep("deepk  n<=12 k<=6 w<=8", 4000, 0x41A7u, 12, 6, 20, 1, 8);
    sweep("deepk  n<=12 k<=6 dense", 4000, 0x41A8u, 12, 6, 60, 1, 8);

    return cpf::report();
}
