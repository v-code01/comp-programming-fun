#pragma once

#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

namespace p0041 {

// undirected weighted edge. weights are int64 because the DP adds them up and
// the sum of a 99-edge tree at w = 1e6 is 9.9e7 -- comfortable, but the
// intermediate INF arithmetic below is what actually forces the width.
struct Edge {
    int u;
    int v;
    std::int64_t w;
};

// "unreachable". picked at 4e17 for one reason: the merge transition adds two dp
// values, and kInf + kInf = 8e17 stays inside int64 (9.22e18). so no guard, no
// branch, no overflow -- and every dp cell is initialized to kInf and only ever
// min'd downward, so nothing can climb past it. real answers top out near 1e8,
// nine orders below, so kInf can never be mistaken for a tree.
constexpr std::int64_t kInf = 400000000000000000LL;

// MINIMUM STEINER TREE -- subset DP over terminals, Dijkstra for the grow step.
//
// vertices are 0-indexed. `terminals` is the set that must end up connected;
// any other vertex may be used as an intermediate (a Steiner point) and often
// must be. returns the minimum total edge weight. k <= 1 -> 0, nothing to join.
//
// -- WHY EXPONENTIAL AT ALL --
// Steiner tree in graphs is NP-hard (Karp '72; MST is the k = n case and is
// easy, shortest path is the k = 2 case and is easy, everything between is not).
// so nobody gets a polynomial algorithm in k. what Dreyfus-Wagner buys is that
// the exponential is confined to k ALONE: O(3^k * n + 2^k * m log n). polynomial
// in the graph, exponential only in the terminal count. that separation is the
// whole point of the problem, and it's why k <= 10 sits in the constraints while
// n runs to 100.
//
// -- LOWER BOUND, HONESTLY --
// you must at least read the graph: Omega(n + m). that's all anyone can prove
// unconditionally. no unconditional superlinear lower bound is known for this
// (nobody can prove one for any problem in NP), so the honest statement is:
// linear is the floor, NP-hardness says you won't hit it, and 3^k / 2^k in the
// terminal count is the best anyone knows -- with a matching conditional lower
// bound under SETH that you can't do 2^{o(k)} * poly(n).
//
// -- THE STATE --
// dp[mask][v] = min weight of a tree that spans the terminals in `mask` AND
// contains vertex v. v need not be a terminal. base: dp[1<<i][terminal_i] = 0,
// the one-vertex tree.
//
// -- THE TWO TRANSITIONS, AND WHY ONE PASS PER MASK IS ENOUGH --
// take an optimal tree T for (mask, v). look at v's degree inside T:
//
//   deg(v) >= 2:  cut T at v. the components hanging off v split into two
//                 nonempty halves; each half plus v is a tree containing v, and
//                 their terminal sets partition mask into sub and mask^sub, both
//                 nonempty. so
//                     dp[mask][v] = dp[sub][v] + dp[mask^sub][v].
//                 both arguments read STRICTLY SMALLER masks -- already final.
//                 that's the MERGE. O(3^k * n) over all masks.
//
//   deg(v) == 1:  v hangs off exactly one edge (v,u). delete it and you have a
//                 tree for the same mask containing u:
//                     dp[mask][v] = dp[mask][u] + w(u,v).
//                 same mask on both sides -- this is a fixpoint, not a scan. it
//                 chains: a Steiner point three edges from the nearest terminal
//                 gets its value through three grow steps. with nonnegative
//                 weights that fixpoint IS a shortest-path closure, so run one
//                 DIJKSTRA over the whole vertex set seeded with every finite
//                 dp[mask][*]. that's the GROW. O(2^k * m log n).
//
//   deg(v) == 0:  T is the single vertex v, so mask is one terminal and that
//                 terminal is v -- the base case.
//
// so: masks in increasing order; inside a mask, MERGE first (it only consumes
// smaller, finished masks), then DIJKSTRA (it closes the same-mask chain). the
// order is a valid topological order of the dependency graph, and it's the
// bug-prone spot: grow before merge and the freshly merged values never get to
// propagate along edges -- you'd under-count every tree whose branch point sits
// away from where the merge happened.
//
// answer = min over v of dp[full][v]. the optimum is a tree (drop any cycle edge
// and it stays connected and gets lighter -- weights are positive), and every
// tree has a vertex, so some v carries it.
//
// complexity O(3^k * n + 2^k * m log n), memory O(2^k * n).
inline std::int64_t steiner_tree(int n, const std::vector<Edge>& edges,
                                 const std::vector<int>& terminals) {
    const int k = static_cast<int>(terminals.size());
    if (k <= 1) return 0;  // zero or one terminal is already connected. 0.

    // adjacency, built once. parallel edges and self-loops are harmless here --
    // a self-loop relaxes v against itself with w >= 0 and never wins.
    std::vector<std::vector<std::pair<int, std::int64_t>>> adj(
        static_cast<std::size_t>(n));
    for (const Edge& e : edges) {
        adj[static_cast<std::size_t>(e.u)].emplace_back(e.v, e.w);
        adj[static_cast<std::size_t>(e.v)].emplace_back(e.u, e.w);
    }

    const int full = (1 << k) - 1;
    // flat [mask * n + v]. the merge's inner loop runs over v with mask fixed, so
    // v is the contiguous axis -- that loop is a straight vectorizable min over
    // three arrays.
    std::vector<std::int64_t> dp(static_cast<std::size_t>(full + 1) *
                                     static_cast<std::size_t>(n),
                                 kInf);
    for (int i = 0; i < k; ++i) {
        dp[static_cast<std::size_t>(1 << i) * static_cast<std::size_t>(n) +
           static_cast<std::size_t>(terminals[static_cast<std::size_t>(i)])] = 0;
    }

    using Item = std::pair<std::int64_t, int>;  // (dist, vertex)
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    for (int mask = 1; mask <= full; ++mask) {
        std::int64_t* row = &dp[static_cast<std::size_t>(mask) *
                                static_cast<std::size_t>(n)];

        // ---- MERGE: two subtrees meeting at v ----
        // submasks of mask, skipping mask itself and the empty set. the pair
        // (sub, mask^sub) is symmetric, so take each pair once: sub < mask^sub.
        // that halves the 3^k work and loses nothing.
        for (int sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask) {
            const int other = mask ^ sub;
            if (sub > other) continue;  // same pair, other order. already did it.
            const std::int64_t* a = &dp[static_cast<std::size_t>(sub) *
                                        static_cast<std::size_t>(n)];
            const std::int64_t* b = &dp[static_cast<std::size_t>(other) *
                                        static_cast<std::size_t>(n)];
            for (int v = 0; v < n; ++v) {
                // kInf + kInf fits, so no reachability branch. a cell that stays
                // at kInf just never wins the min.
                const std::int64_t s = a[v] + b[v];
                if (s < row[v]) row[v] = s;
            }
        }

        // ---- GROW: relax along edges, mask fixed ----
        // multi-source Dijkstra seeded with every finite cell of the row. weights
        // are positive, so the popped-min-is-final invariant holds and one pass
        // settles the whole same-mask fixpoint. lazy deletion instead of a
        // decrease-key -- stale entries are cheap to drop and the heap stays a
        // plain binary heap.
        for (int v = 0; v < n; ++v) {
            if (row[v] < kInf) pq.emplace(row[v], v);
        }
        while (!pq.empty()) {
            const auto [d, v] = pq.top();
            pq.pop();
            if (d > row[v]) continue;  // stale copy -- a better one already ran.
            for (const auto& [u, w] : adj[static_cast<std::size_t>(v)]) {
                const std::int64_t nd = d + w;
                if (nd < row[u]) {
                    row[u] = nd;
                    pq.emplace(nd, u);
                }
            }
        }
    }

    const std::int64_t* last = &dp[static_cast<std::size_t>(full) *
                                   static_cast<std::size_t>(n)];
    std::int64_t best = kInf;
    for (int v = 0; v < n; ++v) {
        if (last[v] < best) best = last[v];
    }
    return best;  // graph is connected, so this is always a real tree weight.
}

}  // namespace p0041
