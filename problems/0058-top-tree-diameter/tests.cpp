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

using p0058::Diam;
using p0058::DynamicDiameter;
using p0058::Reference;
using p0058::TopForest;

// ===========================================================================
// PART 1 -- the top-tree machinery, tested ALONE. before the diameter merge goes
// on, the splay/expose/link/cut/rake plumbing has to be trustworthy. so drive it
// with the SIMPLEST cluster that still exercises compress and rake -- a path length
// -- and gate it on a brute path-walk oracle. if this part is wrong, no merge on
// top of it can be right.
// ===========================================================================

// plain path-length cluster: fold.len is the summed edge weight between the two
// boundaries; rake drops the hanging branch, so a query reads the path, not the
// tree. this isolates the STRUCTURE from the diameter arithmetic.
struct PathLen {
    struct F {
        std::int64_t len = 0;
    };
    static F id() { return F{0}; }
    static F base(std::int64_t w) { return F{w}; }
    static F comp(const F& a, const F& b) { return F{a.len + b.len}; }
    static F rake(const F& a, const F& b) {
        (void)b;
        return F{a.len};
    }
    static void rev(F&) {}
};

// brute path oracle: adjacency + BFS the unique u..v path, sum its weights.
struct BrutePath {
    int n;
    std::vector<std::vector<std::pair<int, std::int64_t>>> adj;
    explicit BrutePath(int k) : n(k), adj(static_cast<std::size_t>(k + 1)) {}
    void link(int u, int v, std::int64_t w) {
        adj[static_cast<std::size_t>(u)].push_back({v, w});
        adj[static_cast<std::size_t>(v)].push_back({u, w});
    }
    void cut(int u, int v) {
        rm(u, v);
        rm(v, u);
    }
    void rm(int a, int b) {
        auto& L = adj[static_cast<std::size_t>(a)];
        for (std::size_t i = 0; i < L.size(); ++i)
            if (L[i].first == b) {
                L.erase(L.begin() + static_cast<std::ptrdiff_t>(i));
                return;
            }
    }
    std::int64_t path(int u, int v) {
        std::vector<int> par(static_cast<std::size_t>(n + 1), 0);
        std::vector<std::int64_t> pw(static_cast<std::size_t>(n + 1), 0);
        std::vector<char> seen(static_cast<std::size_t>(n + 1), 0);
        std::vector<int> q{u};
        seen[static_cast<std::size_t>(u)] = 1;
        for (std::size_t i = 0; i < q.size(); ++i)
            for (auto [w, wt] : adj[static_cast<std::size_t>(q[i])]) {
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                par[static_cast<std::size_t>(w)] = q[i];
                pw[static_cast<std::size_t>(w)] = wt;
                q.push_back(w);
            }
        std::int64_t s = 0;
        for (int c = v; c != u; c = par[static_cast<std::size_t>(c)])
            s += pw[static_cast<std::size_t>(c)];
        return s;
    }
};

// small adjacency mirror the generators read to stay legal: links join different
// trees, cuts hit live edges, queries touch connected pairs.
struct Mirror {
    int n;
    std::vector<std::vector<int>> adj;
    explicit Mirror(int k) : n(k), adj(static_cast<std::size_t>(k + 1)) {}
    std::vector<int> comp(int u) const {
        std::vector<char> seen(static_cast<std::size_t>(n + 1), 0);
        std::vector<int> q{u};
        seen[static_cast<std::size_t>(u)] = 1;
        for (std::size_t i = 0; i < q.size(); ++i)
            for (int w : adj[static_cast<std::size_t>(q[i])])
                if (!seen[static_cast<std::size_t>(w)]) {
                    seen[static_cast<std::size_t>(w)] = 1;
                    q.push_back(w);
                }
        return q;
    }
    bool same(int a, int b) const {
        if (a == b) return true;
        for (int c : comp(a))
            if (c == b) return true;
        return false;
    }
    void link(int a, int b) {
        adj[static_cast<std::size_t>(a)].push_back(b);
        adj[static_cast<std::size_t>(b)].push_back(a);
    }
    void cut(int a, int b) {
        auto rm = [&](int x, int y) {
            auto& L = adj[static_cast<std::size_t>(x)];
            for (std::size_t i = 0; i < L.size(); ++i)
                if (L[i] == y) {
                    L.erase(L.begin() + static_cast<std::ptrdiff_t>(i));
                    break;
                }
        };
        rm(a, b);
        rm(b, a);
    }
    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> e;
        for (int a = 1; a <= n; ++a)
            for (int b : adj[static_cast<std::size_t>(a)])
                if (a < b) e.emplace_back(a, b);
        return e;
    }
};

// one valid random program against the machinery: legal links (across trees, random
// weight), legal cuts (live edges), and path-length queries on connected distinct
// pairs. cuts fire constantly so expose and the rake splices get hammered -- the
// exact spot the structure bugs hide. returns false and prints on the first diff.
bool machinery_one(std::uint64_t seed) {
    cpf::Rng rng(seed);
    int n = static_cast<int>(rng.in_range(1, 12));
    TopForest<PathLen> T(n);
    BrutePath B(n);
    Mirror M(n);
    int m = static_cast<int>(rng.in_range(1, 60));
    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 35) {
            for (int tries = 0; tries < 40; ++tries) {
                int a = static_cast<int>(rng.in_range(1, n));
                int b = static_cast<int>(rng.in_range(1, n));
                if (a != b && !M.same(a, b)) {
                    std::int64_t w = rng.in_range(1, 1000000000LL);
                    M.link(a, b);
                    T.link(a, b, w);
                    B.link(a, b, w);
                    break;
                }
            }
        } else if (roll < 55) {
            auto es = M.edges();
            if (!es.empty()) {
                auto e = es[static_cast<std::size_t>(
                    rng.in_range(0, static_cast<std::int64_t>(es.size()) - 1))];
                M.cut(e.first, e.second);
                T.cut(e.first, e.second);
                B.cut(e.first, e.second);
            }
        } else {
            int u = static_cast<int>(rng.in_range(1, n));
            std::vector<int> c = M.comp(u);
            if (c.size() < 2) continue;  // path needs two distinct endpoints.
            int v = u;
            while (v == u)
                v = c[static_cast<std::size_t>(
                    rng.in_range(0, static_cast<std::int64_t>(c.size()) - 1))];
            if (T.pathAgg(u, v).len != B.path(u, v)) return false;
        }
    }
    return true;
}

// ===========================================================================
// PART 2 -- the diameter top tree vs the BFS/DFS brute oracle.
// ===========================================================================

// one op in the statement's alphabet, vertices 1-based.
//   type 1: link u v w
//   type 2: cut u v
//   type 3: diameter u -> collected
struct Op {
    int type;
    int u;
    int v;
    std::int64_t w;
};

struct Input {
    int n;
    std::vector<Op> ops;
};

// replay a program on anything with the (link, cut, diameter) API, collecting the
// type-3 answers in order. templated so tree and oracle walk the identical program.
template <typename DS>
std::vector<std::int64_t> run(const Input& in) {
    DS ds(in.n);
    std::vector<std::int64_t> out;
    for (const Op& o : in.ops) {
        if (o.type == 1)
            ds.link(o.u, o.v, o.w);
        else if (o.type == 2)
            ds.cut(o.u, o.v);
        else
            out.push_back(ds.diameter(o.u));
    }
    return out;
}

void expect_agree(const Input& in) {
    CPF_CHECK(run<DynamicDiameter>(in) == run<Reference>(in));
}

// a validity-tracking generator. the tree trusts its input, so it keeps an adjacency
// mirror and reads it before every choice: links straddle two trees, cuts hit live
// edges, queries land on any vertex. links and cuts are both biased heavy so the
// tree's SHAPE churns -- paths, stars, caterpillars -- and the diameter keeps moving.
Input gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(1, 14));
    Input in;
    in.n = n;
    Mirror M(n);
    int m = static_cast<int>(rng.in_range(1, 90));
    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 45) {
            for (int tries = 0; tries < 40; ++tries) {
                int a = static_cast<int>(rng.in_range(1, n));
                int b = static_cast<int>(rng.in_range(1, n));
                if (a != b && !M.same(a, b)) {
                    std::int64_t w = rng.in_range(1, 1000000000LL);
                    M.link(a, b);
                    in.ops.push_back({1, a, b, w});
                    break;
                }
            }
        } else if (roll < 65) {
            auto es = M.edges();
            if (!es.empty()) {
                auto e = es[static_cast<std::size_t>(
                    rng.in_range(0, static_cast<std::int64_t>(es.size()) - 1))];
                M.cut(e.first, e.second);
                in.ops.push_back({2, e.first, e.second, 0});
            }
        } else {
            int u = static_cast<int>(rng.in_range(1, n));
            in.ops.push_back({3, u, 0, 0});
        }
    }
    return in;
}

}  // namespace

int main() {
    // ---- PART 1: machinery in isolation, thousands of programs ---------------
    {
        constexpr int kSeqs = 8000;
        int ran = 0, bad = -1;
        for (int s = 0; s < kSeqs; ++s) {
            ++ran;
            if (!machinery_one(0x70F7EE + static_cast<std::uint64_t>(s))) {
                bad = s;
                break;
            }
        }
        CPF_CHECK(bad < 0);
        if (bad >= 0) std::printf("  machinery diff at seq %d\n", bad);
        std::printf("machinery isolation: %d path-length programs vs brute walk\n",
                    ran);
    }

    // ---- hand example: a single edge -> diameter = w -------------------------
    {
        Input in{2, {{1, 1, 2, 7}, {3, 1, 0, 0}, {3, 2, 0, 0}}};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{7, 7}));
        expect_agree(in);
    }

    // ---- hand example: a path -> total weight --------------------------------
    // 1-2-3-4 with weights 3,4,5; diameter is the whole chain 3+4+5 = 12.
    {
        Input in{4,
                 {
                     {1, 1, 2, 3},
                     {1, 2, 3, 4},
                     {1, 3, 4, 5},
                     {3, 1, 0, 0},  // 12
                     {3, 4, 0, 0},  // 12 (same tree)
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{12, 12}));
        expect_agree(in);
    }

    // ---- hand example: a star -> the two heaviest spokes ---------------------
    // center 1, spokes to 2,3,4,5 of weight 10,2,8,5; diameter = 10 + 8 = 18.
    {
        Input in{5,
                 {
                     {1, 1, 2, 10},
                     {1, 1, 3, 2},
                     {1, 1, 4, 8},
                     {1, 1, 5, 5},
                     {3, 3, 0, 0},  // 18
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{18}));
        expect_agree(in);
    }

    // ---- hand example: link two trees, diameter crosses the new edge ---------
    // path 1-2-3 (weights 5,5) has diameter 10; path 4-5-6 (weights 5,5) too.
    // bridge 3-4 with weight 100: the new diameter runs 1..6 = 5+5+100+5+5 = 120.
    {
        Input in{6,
                 {
                     {1, 1, 2, 5}, {1, 2, 3, 5},  // left path, diam 10
                     {1, 4, 5, 5}, {1, 5, 6, 5},  // right path, diam 10
                     {3, 1, 0, 0},                // 10
                     {3, 6, 0, 0},                // 10
                     {1, 3, 4, 100},              // bridge
                     {3, 1, 0, 0},                // 120, crosses the bridge
                     {3, 5, 0, 0},                // 120 (same tree)
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) ==
                   std::vector<std::int64_t>{10, 10, 120, 120}));
        expect_agree(in);
    }

    // ---- hand example: cut, then re-query both halves ------------------------
    // build the 120-diameter tree above, cut the bridge, each half falls to 10.
    {
        Input in{6,
                 {
                     {1, 1, 2, 5}, {1, 2, 3, 5},
                     {1, 4, 5, 5}, {1, 5, 6, 5},
                     {1, 3, 4, 100},
                     {3, 1, 0, 0},  // 120
                     {2, 3, 4, 0},  // cut the bridge
                     {3, 1, 0, 0},  // 10
                     {3, 4, 0, 0},  // 10
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) ==
                   std::vector<std::int64_t>{120, 10, 10}));
        expect_agree(in);
    }

    // ---- edge: a single vertex -> 0 ------------------------------------------
    {
        Input in{1, {{3, 1, 0, 0}}};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{0}));
        expect_agree(in);
    }

    // ---- edge: link then immediate cut -> back to two lone vertices ----------
    {
        Input in{2,
                 {
                     {1, 1, 2, 9},
                     {3, 1, 0, 0},  // 9
                     {2, 1, 2, 0},  // cut
                     {3, 1, 0, 0},  // 0
                     {3, 2, 0, 0},  // 0
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{9, 0, 0}));
        expect_agree(in);
    }

    // ---- edge: a caterpillar -------------------------------------------------
    // spine 1-2-3-4 (weights 2,2,2), legs 5,6,7 hanging off 2,3,4 (weights 9,1,9).
    // diameter runs leg5 -> up to 2 -> along spine to 4 -> down leg7:
    //   9 + 2 + 2 + 9 = 22.
    {
        Input in{7,
                 {
                     {1, 1, 2, 2}, {1, 2, 3, 2}, {1, 3, 4, 2},  // spine
                     {1, 2, 5, 9}, {1, 3, 6, 1}, {1, 4, 7, 9},  // legs
                     {3, 1, 0, 0},                              // 22
                 }};
        CPF_CHECK((run<DynamicDiameter>(in) == std::vector<std::int64_t>{22}));
        expect_agree(in);
    }

    // ---- edge: a big chain, diameter is its full weight ----------------------
    {
        constexpr int kChain = 300;
        Input in;
        in.n = kChain;
        std::int64_t total = 0;
        for (int i = 1; i < kChain; ++i) {
            in.ops.push_back({1, i, i + 1, i});  // edge i has weight i.
            total += i;
        }
        in.ops.push_back({3, 1, 0, 0});
        in.ops.push_back({3, kChain, 0, 0});
        std::vector<std::int64_t> got = run<DynamicDiameter>(in);
        CPF_CHECK((got == std::vector<std::int64_t>{total, total}));
        expect_agree(in);
    }

    // ---- randomized differential vs the brute oracle -------------------------
    // thousands of valid programs. heavy links and cuts keep the tree's shape moving
    // -- paths, stars, caterpillars, forests of many trees -- so the diameter's
    // endpoints wander and the compress+rake merges are stressed from every angle.
    // every type-3 answer is compared; one disagreement prints the seed that replays.
    {
        constexpr int kSeqs = 12000;
        bool ok = cpf::differential(
            kSeqs, 0x0058D1A,
            [](cpf::Rng& rng) { return gen(rng); },
            [](const Input& in) { return run<DynamicDiameter>(in); },
            [](const Input& in) { return run<Reference>(in); });
        CPF_CHECK(ok);
        std::printf("diameter differential: %d valid programs vs brute oracle\n",
                    kSeqs);
    }

    return cpf::report();
}
