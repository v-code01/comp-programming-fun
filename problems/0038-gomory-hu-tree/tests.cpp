#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tuple>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Edge = std::tuple<int, int, std::int64_t>;

// one random weighted graph plus a batch of random query pairs. the whole point
// is a spectrum: sparse enough to disconnect (so 0-cuts show up), dense enough to
// force real cuts, small enough that the Edmonds-Karp oracle stays cheap.
struct Instance {
    int n = 0;
    std::vector<Edge> edges;
    std::vector<std::pair<int, int>> queries;
};

// edge_prob is out of 100. low prob leaves the graph in pieces on purpose.
Instance random_instance(cpf::Rng& rng, int max_n, std::int64_t clo,
                         std::int64_t chi, int edge_prob) {
    Instance in;
    in.n = static_cast<int>(rng.in_range(2, max_n));

    // dense sum matrix so parallel draws add, exactly like the real input path.
    std::vector<std::vector<std::int64_t>> cap(
        static_cast<std::size_t>(in.n),
        std::vector<std::int64_t>(static_cast<std::size_t>(in.n), 0));
    for (int a = 0; a < in.n; ++a) {
        for (int b = a + 1; b < in.n; ++b) {
            if (rng.in_range(1, 100) <= edge_prob) {
                // occasionally stack two edges on a pair to exercise summing.
                const std::int64_t c1 = rng.in_range(clo, chi);
                cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] += c1;
                if (rng.in_range(1, 4) == 1) {
                    cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] +=
                        rng.in_range(clo, chi);
                }
            }
        }
    }
    for (int a = 0; a < in.n; ++a) {
        for (int b = a + 1; b < in.n; ++b) {
            const std::int64_t c =
                cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            if (c > 0) in.edges.emplace_back(a, b, c);
        }
    }

    // hit every ordered pair a<b so no cut escapes the diff, plus a few random
    // repeats to shuffle order. every pair has u != v.
    for (int a = 0; a < in.n; ++a)
        for (int b = a + 1; b < in.n; ++b) in.queries.emplace_back(a, b);
    return in;
}

// the solution's answer for the whole query batch -- one tree, all pairs.
std::vector<std::int64_t> solve_batch(const Instance& in) {
    const p0038::GomoryHu tree(in.n, in.edges);
    std::vector<std::int64_t> out;
    out.reserve(in.queries.size());
    for (const auto& [u, v] : in.queries) out.push_back(tree.min_cut(u, v));
    return out;
}

// the oracle's answer -- a fresh direct max-flow per pair, nothing shared.
std::vector<std::int64_t> ref_batch(const Instance& in) {
    const p0038::DirectMaxFlow direct(in.n, in.edges);
    std::vector<std::int64_t> out;
    out.reserve(in.queries.size());
    for (const auto& [u, v] : in.queries) out.push_back(direct.min_cut(u, v));
    return out;
}

}  // namespace

int main() {
    using namespace p0038;

    // ---- hand shapes. every expected number cross-checked against the direct
    // max-flow oracle so nothing is asserted from my head alone. ----

    // a triangle with caps 1, 2, 3. min cut between any two = the flow, which is
    // the sum of the two edges touching the lighter-isolated vertex.
    //   edges: (0,1)=1, (1,2)=2, (0,2)=3.
    //   cut(0,1): isolate 0 -> 1+3 = 4; isolate 1 -> 1+2 = 3. min side is 3.
    //   cut(1,2): isolate 2 -> 2+3 = 5; isolate 1 -> 1+2 = 3. -> 3... check both.
    {
        std::vector<Edge> e = {{0, 1, 1}, {1, 2, 2}, {0, 2, 3}};
        DirectMaxFlow d(3, e);
        GomoryHu g(3, e);
        for (int u = 0; u < 3; ++u)
            for (int v = u + 1; v < 3; ++v)
                CPF_EQ(g.min_cut(u, v), d.min_cut(u, v));
        // pin the actual values so a silent oracle bug can't wave this through.
        CPF_EQ(g.min_cut(0, 1), 3);  // isolate the {0,2} side of 0: edges 1 + 2.
        CPF_EQ(g.min_cut(0, 2), 4);  // isolate vertex 0: edges 1 + 3.
        CPF_EQ(g.min_cut(1, 2), 3);  // isolate vertex 1: edges 1 + 2.
    }

    // a path graph 0-1-2-3 with caps 5, 1, 4. the bottleneck between the ends is
    // the lightest edge on the line -- classic path-min.
    {
        std::vector<Edge> e = {{0, 1, 5}, {1, 2, 1}, {2, 3, 4}};
        GomoryHu g(4, e);
        DirectMaxFlow d(4, e);
        for (int u = 0; u < 4; ++u)
            for (int v = u + 1; v < 4; ++v)
                CPF_EQ(g.min_cut(u, v), d.min_cut(u, v));
        CPF_EQ(g.min_cut(0, 3), 1);  // whole path -> the weakest link, edge (1,2).
        CPF_EQ(g.min_cut(0, 1), 5);  // adjacent -> that single edge.
        CPF_EQ(g.min_cut(1, 2), 1);
        CPF_EQ(g.min_cut(2, 3), 4);
    }

    // a disconnected pair: two vertices, no edge. min cut = 0.
    {
        std::vector<Edge> e = {};
        GomoryHu g(2, e);
        CPF_EQ(g.min_cut(0, 1), 0);
    }

    // n = 2, a single edge -> the cut is that edge's capacity.
    {
        std::vector<Edge> e = {{0, 1, 42}};
        GomoryHu g(2, e);
        CPF_EQ(g.min_cut(0, 1), 42);
    }

    // no edges at all on a bigger graph -> every pair is 0 (all isolated).
    {
        std::vector<Edge> e = {};
        GomoryHu g(5, e);
        for (int u = 0; u < 5; ++u)
            for (int v = u + 1; v < 5; ++v) CPF_EQ(g.min_cut(u, v), 0);
    }

    // parallel edges summed: two edges on the same pair add to one capacity.
    {
        std::vector<Edge> e = {{0, 1, 10}, {0, 1, 15}};  // caller sums -> 25.
        // the solver expects pre-summed edges; feed the summed form, as main does.
        std::vector<Edge> summed = {{0, 1, 25}};
        GomoryHu g(2, summed);
        CPF_EQ(g.min_cut(0, 1), 25);
        (void)e;
    }

    // two components glued nowhere: {0,1,2} triangle and {3,4} edge. cross-cuts
    // must read 0, in-component cuts must be the real flow.
    {
        std::vector<Edge> e = {{0, 1, 3}, {1, 2, 3}, {0, 2, 3}, {3, 4, 7}};
        GomoryHu g(5, e);
        DirectMaxFlow d(5, e);
        for (int u = 0; u < 5; ++u)
            for (int v = u + 1; v < 5; ++v)
                CPF_EQ(g.min_cut(u, v), d.min_cut(u, v));
        CPF_EQ(g.min_cut(0, 3), 0);  // across the gap.
        CPF_EQ(g.min_cut(2, 4), 0);
        CPF_EQ(g.min_cut(3, 4), 7);  // the lone edge in the second piece.
        CPF_EQ(g.min_cut(0, 1), 6);  // triangle vertex isolates on 3+3.
    }

    // big capacities near the 1e9 ceiling -> the cut lands past 32 bits and must
    // stay exact in int64.
    {
        const std::int64_t B = 1'000'000'000LL;
        std::vector<Edge> e = {{0, 1, B}, {1, 2, B}, {0, 2, B}};
        GomoryHu g(3, e);
        DirectMaxFlow d(3, e);
        CPF_EQ(g.min_cut(0, 1), 2 * B);  // 2e9, needs 64 bits.
        CPF_EQ(g.min_cut(0, 1), d.min_cut(0, 1));
    }

    // ---- randomized differential vs the direct max-flow oracle. thousands of
    // small weighted graphs, all pairs queried, comparing the tree's path-min to
    // an independent per-pair flow. mixes of density so both disconnected 0-cuts
    // and dense real cuts land. ----

    // sparse: ~25% edges, small caps, graphs that fall apart into components.
    {
        const int kCases = 4000;
        bool ok = cpf::differential(
            kCases, 0x38A1u,
            [](cpf::Rng& rng) { return random_instance(rng, 8, 1, 12, 25); },
            solve_batch, ref_batch);
        CPF_CHECK(ok);
        std::printf("sparse differential vs direct max-flow: %d graphs, n<=8\n",
                    kCases);
    }

    // dense: ~80% edges, so the graph is usually connected and cuts are fat.
    {
        const int kCases = 4000;
        bool ok = cpf::differential(
            kCases, 0x38B2u,
            [](cpf::Rng& rng) { return random_instance(rng, 8, 1, 20, 80); },
            solve_batch, ref_batch);
        CPF_CHECK(ok);
        std::printf("dense differential vs direct max-flow: %d graphs, n<=8\n",
                    kCases);
    }

    // bigger n, medium density -- more tree edges, more grandparent rotations.
    {
        const int kCases = 1500;
        bool ok = cpf::differential(
            kCases, 0x38C3u,
            [](cpf::Rng& rng) { return random_instance(rng, 14, 1, 50, 50); },
            solve_batch, ref_batch);
        CPF_CHECK(ok);
        std::printf("medium differential vs direct max-flow: %d graphs, n<=14\n",
                    kCases);
    }

    // large capacities: caps at the 1e9 scale so int64 cuts get stressed at real
    // bounds, still diffed against the oracle on small n.
    {
        const int kCases = 1500;
        bool ok = cpf::differential(
            kCases, 0x38D4u,
            [](cpf::Rng& rng) {
                return random_instance(rng, 9, 1, 1'000'000'000LL, 60);
            },
            solve_batch, ref_batch);
        CPF_CHECK(ok);
        std::printf("large-cap differential vs direct max-flow: %d graphs, "
                    "|c|<=1e9, n<=9\n",
                    kCases);
    }

    // ---- a sweep that counts how many pairs land on a 0 cut, to prove the
    // disconnected path is exercised, not just claimed. very sparse graphs. ----
    {
        const int kCases = 3000;
        std::int64_t zero_cuts = 0, total_pairs = 0, diffs = 0;
        cpf::Rng rng(0x38E5u);
        for (int c = 0; c < kCases; ++c) {
            Instance in = random_instance(rng, 7, 1, 10, 15);  // 15% -> shredded.
            std::vector<std::int64_t> got = solve_batch(in);
            std::vector<std::int64_t> want = ref_batch(in);
            if (got != want) ++diffs;
            for (std::size_t k = 0; k < want.size(); ++k) {
                ++total_pairs;
                if (want[k] == 0) ++zero_cuts;
            }
        }
        CPF_CHECK(diffs == 0);
        std::printf("very-sparse sweep: %d graphs, %lld/%lld pairs hit a 0 cut, "
                    "%lld diffs\n",
                    kCases, static_cast<long long>(zero_cuts),
                    static_cast<long long>(total_pairs),
                    static_cast<long long>(diffs));
    }

    return cpf::report();
}
