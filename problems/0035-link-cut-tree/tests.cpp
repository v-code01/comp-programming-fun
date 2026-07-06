#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0035::LinkCutTree;
using p0035::Reference;

// one op in the statement's alphabet, vertices 1-based.
//   type 1: link u v          (x unused)
//   type 2: cut u v           (x unused)
//   type 3: pathAdd u v x
//   type 4: pathSum u v -> collected  (x unused)
struct Op {
    int type;
    int u;
    int v;
    std::int64_t x;
};

struct Input {
    std::vector<std::int64_t> init;
    std::vector<Op> ops;
};

// replay a program on any structure with the LCT API, collecting the type-4
// answers in order. templated so the tree and the oracle walk the identical path.
template <typename DS>
std::vector<std::int64_t> run(DS& ds, const Input& in) {
    ds.build(in.init);
    std::vector<std::int64_t> out;
    for (const Op& o : in.ops) {
        if (o.type == 1)
            ds.link(o.u, o.v);
        else if (o.type == 2)
            ds.cut(o.u, o.v);
        else if (o.type == 3)
            ds.pathAdd(o.u, o.v, o.x);
        else
            out.push_back(ds.pathSum(o.u, o.v));
    }
    return out;
}

// tree vs oracle on one fixed program -- both must emit the same type-4 stream.
void expect_agree(const Input& in) {
    LinkCutTree lct;
    Reference ref;
    CPF_CHECK(run(lct, in) == run(ref, in));
}

// ---- a connectivity-tracking generator -------------------------------------
// the LCT trusts its inputs: links join different trees, cuts hit real edges,
// path ops touch connected pairs. so the generator has to only ever emit legal
// ops -- it keeps its own adjacency mirror and reads it before every choice. n is
// small so the O(n) checks are free and the O(n)-per-op oracle stays cheap.

// the component of u, by BFS over the current adjacency.
std::vector<int> component(const std::vector<std::vector<int>>& adj, int u) {
    std::vector<char> seen(adj.size(), 0);
    std::vector<int> q;
    q.push_back(u);
    seen[static_cast<std::size_t>(u)] = 1;
    for (std::size_t i = 0; i < q.size(); ++i) {
        int x = q[i];
        for (int w : adj[static_cast<std::size_t>(x)]) {
            if (seen[static_cast<std::size_t>(w)]) continue;
            seen[static_cast<std::size_t>(w)] = 1;
            q.push_back(w);
        }
    }
    return q;
}

bool same_comp(const std::vector<std::vector<int>>& adj, int a, int b) {
    if (a == b) return true;
    for (int c : component(adj, a))
        if (c == b) return true;
    return false;
}

Input gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 12));
    Input in;
    in.init.resize(static_cast<std::size_t>(n));
    for (auto& v : in.init) v = rng.in_range(-1000000000LL, 1000000000LL);

    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n + 1));
    int m = static_cast<int>(rng.in_range(1, 40));

    // pick a connected (u, v), u may equal v -- shared by add and query.
    auto connected_pair = [&](int& u, int& v) {
        u = static_cast<int>(rng.in_range(1, n));
        std::vector<int> comp = component(adj, u);
        v = comp[static_cast<std::size_t>(
            rng.in_range(0, static_cast<std::int64_t>(comp.size()) - 1))];
    };

    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));

        if (roll < 30) {
            // link -- sample pairs until one straddles two components.
            for (int tries = 0; tries < 40; ++tries) {
                int a = static_cast<int>(rng.in_range(1, n));
                int b = static_cast<int>(rng.in_range(1, n));
                if (a != b && !same_comp(adj, a, b)) {
                    adj[static_cast<std::size_t>(a)].push_back(b);
                    adj[static_cast<std::size_t>(b)].push_back(a);
                    in.ops.push_back({1, a, b, 0});
                    break;
                }
            }
            continue;  // if no pair was free (one component), skip this slot.
        }

        if (roll < 50) {
            // cut -- gather current edges, drop a random one from both sides.
            std::vector<std::pair<int, int>> edges;
            for (int a = 1; a <= n; ++a)
                for (int b : adj[static_cast<std::size_t>(a)])
                    if (a < b) edges.emplace_back(a, b);
            if (!edges.empty()) {
                auto e = edges[static_cast<std::size_t>(rng.in_range(
                    0, static_cast<std::int64_t>(edges.size()) - 1))];
                auto rm = [&](int a, int b) {
                    auto& L = adj[static_cast<std::size_t>(a)];
                    for (std::size_t i = 0; i < L.size(); ++i)
                        if (L[i] == b) {
                            L.erase(L.begin() + static_cast<std::ptrdiff_t>(i));
                            break;
                        }
                };
                rm(e.first, e.second);
                rm(e.second, e.first);
                in.ops.push_back({2, e.first, e.second, 0});
            }
            continue;
        }

        if (roll < 75) {
            // pathAdd on a connected pair, signed x.
            int u = 0, v = 0;
            connected_pair(u, v);
            std::int64_t x = rng.in_range(-1000000000LL, 1000000000LL);
            in.ops.push_back({3, u, v, x});
            continue;
        }

        // pathSum on a connected pair.
        int u = 0, v = 0;
        connected_pair(u, v);
        in.ops.push_back({4, u, v, 0});
    }

    return in;
}

}  // namespace

int main() {
    // ---- hand example: a chain, add on a subpath, cut, re-query -------------
    // vertices 1..5 with values 1 2 3 4 5, linked into the chain 1-2-3-4-5.
    {
        Input in;
        in.init = {1, 2, 3, 4, 5};
        in.ops = {
            {1, 1, 2, 0},      // link 1-2
            {1, 2, 3, 0},      // link 2-3
            {1, 3, 4, 0},      // link 3-4
            {1, 4, 5, 0},      // link 4-5
            {4, 1, 5, 0},      // sum 1..5 -> 15
            {4, 2, 4, 0},      // sum 2..4 -> 9
            {3, 2, 4, 10},     // add 10 on 2..4 -> 1 12 13 14 5
            {4, 1, 5, 0},      // sum 1..5 -> 45
            {4, 3, 3, 0},      // sum 3..3 -> 13
            {2, 2, 3, 0},      // cut 2-3 -> {1,2} {3,4,5}
            {4, 4, 5, 0},      // sum 4..5 -> 19
            {4, 1, 2, 0},      // sum 1..2 -> 13
            {3, 3, 5, -100},   // add -100 on 3..5 -> -87 -86 -95
            {4, 3, 5, 0},      // sum 3..5 -> -268
        };
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        std::vector<std::int64_t> want{15, 9, 45, 13, 19, 13, -268};
        CPF_CHECK(got == want);
        expect_agree(in);
    }

    // ---- hand example: a star, path crosses the center ---------------------
    // center 1, leaves 2 3 4; values 10 20 30 40.
    {
        Input in;
        in.init = {10, 20, 30, 40};
        in.ops = {
            {1, 1, 2, 0},   // link
            {1, 1, 3, 0},
            {1, 1, 4, 0},
            {4, 2, 3, 0},   // path 2-1-3 -> 20+10+30 = 60
            {4, 4, 4, 0},   // single vertex -> 40
            {3, 2, 4, 5},   // add 5 on 2-1-4 -> 25 15 30 45
            {4, 2, 3, 0},   // 25+15+30 = 70
            {4, 3, 4, 0},   // 30+15+45 = 90
        };
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        std::vector<std::int64_t> want{60, 40, 70, 90};
        CPF_CHECK(got == want);
        expect_agree(in);
    }

    // ---- enumerated edges --------------------------------------------------
    // single vertex, u == v: the path is one node, sum is its value, add hits it.
    {
        Input in;
        in.init = {7};
        in.ops = {
            {4, 1, 1, 0},   // sum -> 7
            {3, 1, 1, 100}, // add 100 -> 107
            {4, 1, 1, 0},   // sum -> 107
        };
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        CPF_CHECK((got == std::vector<std::int64_t>{7, 107}));
        expect_agree(in);
    }

    // link then immediately cut -- back to two isolated vertices.
    {
        Input in;
        in.init = {3, 5};
        in.ops = {
            {1, 1, 2, 0},   // link
            {4, 1, 2, 0},   // 3+5 = 8
            {2, 1, 2, 0},   // cut
            {4, 1, 1, 0},   // 3
            {4, 2, 2, 0},   // 5
        };
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        CPF_CHECK((got == std::vector<std::int64_t>{8, 3, 5}));
        expect_agree(in);
    }

    // a long chain: link 1..200, add on the whole thing, sum end to end, then
    // re-root implicitly by querying the reverse direction.
    {
        constexpr int kChain = 200;
        Input in;
        in.init.assign(static_cast<std::size_t>(kChain), 1);  // all ones.
        for (int i = 1; i < kChain; ++i) in.ops.push_back({1, i, i + 1, 0});
        in.ops.push_back({4, 1, kChain, 0});        // sum -> 200
        in.ops.push_back({4, kChain, 1, 0});        // reversed query -> 200
        in.ops.push_back({3, 1, kChain, 3});        // add 3 to all -> each is 4
        in.ops.push_back({4, 50, 150, 0});          // 101 vertices * 4 = 404
        in.ops.push_back({3, 1, 100, -10});         // vertices 1..100 -> -6
        in.ops.push_back({4, 1, kChain, 0});        // 100*(-6) + 100*4 = -200
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        std::vector<std::int64_t> want{kChain, kChain, 404, -200};
        CPF_CHECK(got == want);
        expect_agree(in);
    }

    // a forest that never fully connects -- two separate chains, queried apart.
    {
        Input in;
        in.init = {1, 2, 3, 4, 5, 6};
        in.ops = {
            {1, 1, 2, 0}, {1, 2, 3, 0},   // chain 1-2-3
            {1, 4, 5, 0}, {1, 5, 6, 0},   // chain 4-5-6
            {4, 1, 3, 0},                 // 1+2+3 = 6
            {4, 4, 6, 0},                 // 4+5+6 = 15
            {3, 4, 6, 1000000000LL},      // add 1e9 to the second chain
            {4, 4, 6, 0},                 // 15 + 3e9 = 3000000015
            {4, 1, 3, 0},                 // first chain untouched -> 6
        };
        LinkCutTree lct;
        std::vector<std::int64_t> got = run(lct, in);
        std::vector<std::int64_t> want{6, 15, 3000000015LL, 6};
        CPF_CHECK(got == want);
        expect_agree(in);
    }

    // ---- randomized differential vs the oracle -----------------------------
    // thousands of independent valid programs. links/cuts churn the preferred-path
    // structure so evert and lazy-reverse fire constantly; interleaved path-adds
    // with signed x stress the reverse+add composition inside push_down; every
    // type-4 answer is compared. one disagreement prints the seed that replays it.
    {
        constexpr int kSeqs = 6000;
        bool ok = cpf::differential(
            kSeqs, 0x11CD7,
            [](cpf::Rng& rng) { return gen(rng); },
            [](const Input& in) {
                LinkCutTree lct;
                return run(lct, in);
            },
            [](const Input& in) {
                Reference ref;
                return run(ref, in);
            });
        CPF_CHECK(ok);
        std::printf("differential: %d valid op-programs vs brute oracle\n", kSeqs);
    }

    return cpf::report();
}
