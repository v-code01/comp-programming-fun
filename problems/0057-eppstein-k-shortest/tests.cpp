#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

// case counts and the oracle pop budget. tunable at the compiler so the sweep
// can be widened without touching the file; the defaults are what ctest runs.
#ifndef CASES1
#define CASES1 25000
#endif
#ifndef CASES2
#define CASES2 25000
#endif
// legitimate tiny cases settle in well under 8k pops; the budget only stops
// zero-cycle starvation, so 40k is a comfortable margin that skips nothing real.
#ifndef ORACLE_CAP
#define ORACLE_CAP 40000
#endif

namespace {

using p0057::Edge;
using Lens = std::vector<std::int64_t>;

// can s reach t at all? forward BFS -- the ground truth for the all-(-1) case,
// which a zero-cycle can stop the best-first oracle from ever settling.
bool reachable(int n, int s, int t, const std::vector<Edge>& edges) {
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n) + 1);
    for (const auto& e : edges) adj[static_cast<std::size_t>(e.u)].push_back(e.v);
    std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> stk{s};
    seen[static_cast<std::size_t>(s)] = 1;
    while (!stk.empty()) {
        const int v = stk.back();
        stk.pop_back();
        if (v == t) return true;
        for (const int x : adj[static_cast<std::size_t>(v)])
            if (!seen[static_cast<std::size_t>(x)]) {
                seen[static_cast<std::size_t>(x)] = 1;
                stk.push_back(x);
            }
    }
    return s == t;  // the empty walk counts when s == t.
}

// a random tiny digraph. edges drawn with replacement, so multi-edges, cycles,
// self-loops, and dead ends all arrive on their own -- every place Eppstein's
// heap-of-heaps could quietly be wrong. weights lean positive with the odd zero,
// so zero-weight cycles (which starve the oracle) stay rare but present.
struct Instance {
    int n, s, t, k;
    std::vector<Edge> edges;
};

Instance random_graph(cpf::Rng& rng, int max_n, int max_m, int max_k,
                      int max_w, int zero_pct) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(1, max_n));
    in.s = static_cast<int>(rng.in_range(1, in.n));
    in.t = static_cast<int>(rng.in_range(1, in.n));
    in.k = static_cast<int>(rng.in_range(1, max_k));
    const int m = static_cast<int>(rng.in_range(0, max_m));
    for (int e = 0; e < m; ++e) {
        const int u = static_cast<int>(rng.in_range(1, in.n));
        const int v = static_cast<int>(rng.in_range(1, in.n));
        const bool zero = rng.in_range(0, 99) < zero_pct;
        const std::int64_t w = zero ? 0 : rng.in_range(0, max_w);
        in.edges.push_back({u, v, w});
    }
    return in;
}

void dump(const Instance& in, const Lens& got, const Lens& want) {
    std::printf("  n=%d s=%d t=%d k=%d edges:", in.n, in.s, in.t, in.k);
    for (const auto& e : in.edges)
        std::printf(" %d->%d(%lld)", e.u, e.v, static_cast<long long>(e.w));
    std::printf("\n  got :");
    for (const auto x : got) std::printf(" %lld", static_cast<long long>(x));
    std::printf("\n  want:");
    for (const auto x : want) std::printf(" %lld", static_cast<long long>(x));
    std::printf("\n");
}

}  // namespace

int main() {
    using namespace p0057;

    // ---- hand shapes. each expected vector is reasoned out AND cross-checked
    // against the naive oracle where the oracle terminates. ----

    // a line 1->2->3->4, one path only. k=3: the single length, then -1, -1.
    {
        std::vector<Edge> e = {{1, 2, 2}, {2, 3, 3}, {3, 4, 5}};
        Lens got = k_shortest_walks(4, 1, 4, 3, e);
        CPF_EQ(got[0], 10);
        CPF_EQ(got[1], -1);
        CPF_EQ(got[2], -1);
    }

    // two parallel paths 1->2->4 (len 2) and 1->3->4 (len 6). k=3 -> both, -1.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 4, 1}, {1, 3, 5}, {3, 4, 1}};
        Lens got = k_shortest_walks(4, 1, 4, 3, e);
        CPF_EQ(got[0], 2);
        CPF_EQ(got[1], 6);
        CPF_EQ(got[2], -1);
    }

    // a positive multi-vertex cycle. 1->3 (len 1); the 1<->2 loop costs 2 a lap,
    // so lengths are 1,3,5,7,... an arithmetic run of +2.
    {
        std::vector<Edge> e = {{1, 2, 1}, {2, 1, 1}, {1, 3, 1}};
        Lens got = k_shortest_walks(3, 1, 3, 4, e);
        CPF_EQ(got[0], 1);
        CPF_EQ(got[1], 3);
        CPF_EQ(got[2], 5);
        CPF_EQ(got[3], 7);
    }

    // t unreachable: 1->2 only, ask for 3. no walk exists -- all -1.
    {
        std::vector<Edge> e = {{1, 2, 4}};
        Lens got = k_shortest_walks(3, 1, 3, 3, e);
        CPF_EQ(got[0], -1);
        CPF_EQ(got[1], -1);
        CPF_EQ(got[2], -1);
    }

    // ---- edges the machinery has to survive. ----

    // s == t. the empty walk is length 0; the 1<->2 loop costs 10 a lap, so
    // 0,10,20.
    {
        std::vector<Edge> e = {{1, 2, 5}, {2, 1, 5}};
        Lens got = k_shortest_walks(2, 1, 1, 3, e);
        CPF_EQ(got[0], 0);
        CPF_EQ(got[1], 10);
        CPF_EQ(got[2], 20);
    }

    // s == t with no way back. only the empty walk, length 0, then -1.
    {
        std::vector<Edge> e = {{1, 2, 5}};
        Lens got = k_shortest_walks(2, 1, 1, 3, e);
        CPF_EQ(got[0], 0);
        CPF_EQ(got[1], -1);
        CPF_EQ(got[2], -1);
    }

    // a self-loop giving an arithmetic progression. 1->2 (len 1), and 1->1 (+3)
    // any number of times first: 1,4,7,10,...
    {
        std::vector<Edge> e = {{1, 1, 3}, {1, 2, 1}};
        Lens got = k_shortest_walks(2, 1, 2, 5, e);
        CPF_EQ(got[0], 1);
        CPF_EQ(got[1], 4);
        CPF_EQ(got[2], 7);
        CPF_EQ(got[3], 10);
        CPF_EQ(got[4], 13);
    }

    // a ZERO-weight self-loop on the path: infinitely many walks all tied. every
    // length equals the one path length. 1->2 (len 4), loop 1->1 at cost 0.
    {
        std::vector<Edge> e = {{1, 1, 0}, {1, 2, 4}};
        Lens got = k_shortest_walks(2, 1, 2, 4, e);
        CPF_EQ(got[0], 4);
        CPF_EQ(got[1], 4);
        CPF_EQ(got[2], 4);
        CPF_EQ(got[3], 4);
    }

    // n = 1, s == t, no edges. exactly the empty walk.
    {
        std::vector<Edge> e = {};
        Lens got = k_shortest_walks(1, 1, 1, 2, e);
        CPF_EQ(got[0], 0);
        CPF_EQ(got[1], -1);
    }

    // parallel edges to t of different weight, so a vertex has a real sidetrack
    // that is NOT its tree edge. 1->2 twice (w 3 and w 7); best walk 3, next 7.
    {
        std::vector<Edge> e = {{1, 2, 3}, {1, 2, 7}};
        Lens got = k_shortest_walks(2, 1, 2, 3, e);
        CPF_EQ(got[0], 3);
        CPF_EQ(got[1], 7);
        CPF_EQ(got[2], -1);
    }

    // a diamond with a detour: two ways into t plus a cycle to spin. exercises a
    // cross-edge landing in a non-empty inherited heap.
    {
        std::vector<Edge> e = {{1, 2, 1}, {1, 3, 1}, {2, 4, 1},
                               {3, 4, 1}, {4, 1, 1}};
        Lens got = k_shortest_walks(4, 1, 4, 3, e);
        ReferenceSolver ref;
        auto r = ref.solve(4, 1, 4, 3, e);
        CPF_CHECK(r.complete);
        CPF_CHECK(got == r.lens);
    }

    // ---- randomized differential vs the naive oracle. thousands of tiny graphs.
    // unreachable cases are checked against forward reachability (ground truth);
    // reachable cases against the oracle when it settles. ----
    {
        const int kCases = CASES1;
        int diffs = 0, skipped = 0, unreachable = 0, with_zero = 0, tied = 0;
        bool first = true;
        cpf::Rng rng(0x0057u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_graph(rng, 7, 12, 10, 6, 12);
            Lens got = k_shortest_walks(in.n, in.s, in.t, in.k, in.edges);

            for (const auto& e : in.edges)
                if (e.w == 0) {
                    ++with_zero;
                    break;
                }

            if (!reachable(in.n, in.s, in.t, in.edges)) {
                ++unreachable;
                bool all_neg = true;
                for (const auto x : got)
                    if (x != -1) all_neg = false;
                if (!all_neg) {
                    ++diffs;
                    if (first) {
                        std::printf("  DIFF (unreachable not all -1) case %d\n", c);
                        Lens want(got.size(), -1);
                        dump(in, got, want);
                        first = false;
                    }
                }
                continue;
            }

            ReferenceSolver ref;
            auto r = ref.solve(in.n, in.s, in.t, in.k, in.edges, ORACLE_CAP);
            if (!r.complete) {
                ++skipped;
                continue;
            }
            for (int i = 0; i + 1 < in.k; ++i)
                if (r.lens[static_cast<std::size_t>(i)] ==
                        r.lens[static_cast<std::size_t>(i + 1)] &&
                    r.lens[static_cast<std::size_t>(i)] != -1) {
                    ++tied;
                    break;
                }
            if (got != r.lens) {
                ++diffs;
                if (first) {
                    std::printf("  DIFF at case %d\n", c);
                    dump(in, got, r.lens);
                    first = false;
                }
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf(
            "differential n<=7 m<=12 k<=10: %d cases, %d diffs, %d skipped "
            "(oracle budget), %d unreachable, %d with a zero edge, %d with tied "
            "lengths\n",
            kCases, diffs, skipped, unreachable, with_zero, tied);
    }

    // ---- a denser, cycle-heavy sweep: more edges on fewer vertices so real
    // sidetracks and thick cycles dominate and the persistent meld gets hammered.
    // weights all >= 1 here, so the oracle rarely starves. ----
    {
        const int kCases = CASES2;
        int diffs = 0, skipped = 0;
        bool first = true;
        cpf::Rng rng(0xD0057u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_graph(rng, 5, 16, 12, 4, 0);
            Lens got = k_shortest_walks(in.n, in.s, in.t, in.k, in.edges);
            if (!reachable(in.n, in.s, in.t, in.edges)) {
                for (const auto x : got) CPF_CHECK(x == -1);
                continue;
            }
            ReferenceSolver ref;
            auto r = ref.solve(in.n, in.s, in.t, in.k, in.edges, ORACLE_CAP);
            if (!r.complete) {
                ++skipped;
                continue;
            }
            if (got != r.lens) {
                ++diffs;
                if (first) {
                    std::printf("  DIFF at dense case %d\n", c);
                    dump(in, got, r.lens);
                    first = false;
                }
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf("dense differential n<=5 m<=16 k<=12: %d cases, %d diffs, "
                    "%d skipped\n",
                    kCases, diffs, skipped);
    }

    return cpf::report();
}
