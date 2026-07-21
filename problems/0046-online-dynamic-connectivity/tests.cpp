#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using p0046::DynamicConnectivity;
using p0046::EulerTourForest;
using p0046::Op;

// ===========================================================================
// PART 1 -- the Euler-tour forest, tested ALONE. before a single level goes on
// top, the link/cut/connected/aggregate layer has to be trustworthy. so gate it
// on its own brute oracle first -- a forest that stores adjacency and answers by
// flood fill. if this part is wrong, nothing above it can be right.
// ===========================================================================

// brute forest: an adjacency mirror. links only ever join two components, cuts
// only ever hit a live edge, marks are a plain per-vertex flag.
struct BruteForest {
    int n = 0;
    std::vector<std::set<int>> adj;
    std::vector<char> mark;

    void init(int k) {
        n = k;
        adj.assign(static_cast<std::size_t>(k) + 1, {});
        mark.assign(static_cast<std::size_t>(k) + 1, 0);
    }
    std::vector<int> comp(int u) const {
        std::vector<char> seen(static_cast<std::size_t>(n) + 1, 0);
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
    bool connected(int u, int v) const {
        if (u == v) return true;
        for (int c : comp(u))
            if (c == v) return true;
        return false;
    }
    int comp_size(int u) const { return static_cast<int>(comp(u).size()); }
    void link(int u, int v) {
        adj[static_cast<std::size_t>(u)].insert(v);
        adj[static_cast<std::size_t>(v)].insert(u);
    }
    void cut(int u, int v) {
        adj[static_cast<std::size_t>(u)].erase(v);
        adj[static_cast<std::size_t>(v)].erase(u);
    }
    bool any_marked(int u) const {
        for (int c : comp(u))
            if (mark[static_cast<std::size_t>(c)]) return true;
        return false;
    }
};

// one valid random program against the ETT: legal links (across components), legal
// cuts (live edges), mark toggles both ways, and queries that check connectivity,
// component size, and the marked-vertex descent all at once. cuts fire constantly,
// so the tour splits and rotations get hammered -- the exact spot the cyclic-rotate
// bug hid. returns false and prints on the first disagreement.
bool ett_one(std::uint64_t seed) {
    cpf::Rng rng(seed);
    int n = static_cast<int>(rng.in_range(1, 12));
    EulerTourForest F;
    F.init(n);
    BruteForest B;
    B.init(n);
    std::map<std::pair<int, int>, EulerTourForest::Edge> handles;
    int edge_id = 1;

    int m = static_cast<int>(rng.in_range(1, 60));
    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 30) {
            for (int tries = 0; tries < 40; ++tries) {
                int a = static_cast<int>(rng.in_range(1, n));
                int b = static_cast<int>(rng.in_range(1, n));
                if (a != b && !B.connected(a, b)) {
                    int lo = a < b ? a : b, hi = a < b ? b : a;
                    handles[{lo, hi}] = F.link(a, b, edge_id++);
                    B.link(a, b);
                    break;
                }
            }
        } else if (roll < 50) {
            if (!handles.empty()) {
                int idx = static_cast<int>(rng.in_range(
                    0, static_cast<std::int64_t>(handles.size()) - 1));
                auto it = handles.begin();
                std::advance(it, idx);
                F.cut(it->second);
                B.cut(it->first.first, it->first.second);
                handles.erase(it);
            }
        } else if (roll < 65) {
            int a = static_cast<int>(rng.in_range(1, n));
            bool on = !B.mark[static_cast<std::size_t>(a)];
            F.set_vertex_mark(a, on);
            B.mark[static_cast<std::size_t>(a)] = on ? 1 : 0;
        } else {
            int a = static_cast<int>(rng.in_range(1, n));
            int b = static_cast<int>(rng.in_range(1, n));
            if (F.connected(a, b) != B.connected(a, b)) return false;
            if (F.comp_size(a) != B.comp_size(a)) return false;
            int mv = F.find_marked_vertex(a);
            if (B.any_marked(a)) {
                // the descent must return SOME marked vertex of a's component.
                if (mv == 0) return false;
                if (!B.connected(a, mv)) return false;
                if (!B.mark[static_cast<std::size_t>(mv)]) return false;
            } else {
                if (mv != 0) return false;  // no mark anywhere -> must be 0.
            }
        }
    }
    return true;
}

// ===========================================================================
// PART 2 -- the full HDLT structure vs the brute recompute oracle.
// ===========================================================================

struct Input {
    int n;
    std::vector<Op> ops;
};

// replay a program on the online structure, collecting the type-3 answers.
std::vector<int> run_hdlt(const Input& in) {
    DynamicConnectivity dc;
    dc.init(in.n);
    std::vector<int> out;
    for (const Op& o : in.ops) {
        if (o.type == 1)
            dc.insert(o.u, o.v);
        else if (o.type == 2)
            dc.remove(o.u, o.v);
        else
            out.push_back(dc.connected(o.u, o.v) ? 1 : 0);
    }
    return out;
}

std::vector<int> run_ref(const Input& in) {
    return p0046::referenceSolve(in.n, in.ops);
}

void expect_agree(const Input& in) { CPF_CHECK(run_hdlt(in) == run_ref(in)); }

// a validity-tracking generator. the structure trusts its input -- inserts hit an
// absent edge, deletes hit a present one -- so the generator carries a live-edge
// mirror and reads it before every choice. inserts are biased high to build dense,
// cyclic graphs; deletes are biased high too, so tree edges get cut and the
// replacement search has non-tree edges to actually find. that mix is the whole
// point: it's what forces the level machinery to run.
Input gen(cpf::Rng& rng) {
    int n = static_cast<int>(rng.in_range(2, 16));
    Input in;
    in.n = n;
    std::set<std::pair<int, int>> present;
    auto norm = [](int u, int v) {
        return u < v ? std::make_pair(u, v) : std::make_pair(v, u);
    };

    int m = static_cast<int>(rng.in_range(1, 140));
    for (int t = 0; t < m; ++t) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 45) {
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
        } else if (roll < 75) {
            if (!present.empty()) {
                int idx = static_cast<int>(rng.in_range(
                    0, static_cast<std::int64_t>(present.size()) - 1));
                auto it = present.begin();
                std::advance(it, idx);
                in.ops.push_back({2, it->first, it->second});
                present.erase(it);
            }
        } else {
            int a = static_cast<int>(rng.in_range(1, n));
            int b = static_cast<int>(rng.in_range(1, n));
            in.ops.push_back({3, a, b});
        }
    }
    return in;
}

}  // namespace

int main() {
    // ---- PART 1: ETT in isolation, thousands of programs -------------------
    {
        constexpr int kSeqs = 4000;
        int ran = 0, bad_seed = -1;
        for (int s = 0; s < kSeqs; ++s) {
            ++ran;
            if (!ett_one(0xE770F + static_cast<std::uint64_t>(s))) {
                bad_seed = s;
                break;
            }
        }
        CPF_CHECK(bad_seed < 0);
        if (bad_seed >= 0) std::printf("  ETT diff at seq %d\n", bad_seed);
        std::printf("ETT isolation differential: %d programs vs brute forest\n",
                    ran);
    }

    // ---- PART 2 hand example: an edge, a query, a delete -------------------
    {
        Input in{2, {{1, 1, 2}, {3, 1, 2}, {2, 1, 2}, {3, 1, 2}}};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 0}));
        expect_agree(in);
    }

    // ---- hand example: THE replacement case --------------------------------
    // a triangle 1-2-3. every pair connected. delete the tree edge that carries the
    // spanning path, and a non-tree edge of the cycle has to be promoted to keep
    // 1 and 3 connected. this is the case the whole level machinery exists for.
    {
        Input in{3,
                 {
                     {1, 1, 2},  // tree edge
                     {1, 2, 3},  // tree edge -- spanning path 1-2-3
                     {1, 1, 3},  // closes the cycle -> non-tree edge
                     {3, 1, 3},  // connected -> 1
                     {2, 2, 3},  // cut a tree edge; 1-3 must survive via replacement
                     {3, 1, 3},  // still connected -> 1
                     {3, 2, 3},  // 2 rejoined through 1 -> 1
                     {2, 1, 3},  // drop the replacement too
                     {3, 2, 3},  // now 2 is cut off from 3 -> 0
                     {3, 1, 2},  // 1-2 still there -> 1
                 }};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 1, 1, 0, 1}));
        expect_agree(in);
    }

    // ---- hand example: two components merge, then split --------------------
    {
        Input in{4,
                 {
                     {1, 1, 2},  // {1,2}
                     {1, 3, 4},  // {3,4}
                     {3, 1, 3},  // apart -> 0
                     {1, 2, 3},  // bridge the two -> one component
                     {3, 1, 4},  // -> 1
                     {2, 2, 3},  // drop the bridge (no replacement exists)
                     {3, 1, 4},  // -> 0
                     {3, 1, 2},  // left half intact -> 1
                     {3, 3, 4},  // right half intact -> 1
                 }};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{0, 1, 0, 1, 1}));
        expect_agree(in);
    }

    // ---- hand example: a chain that breaks in the middle -------------------
    {
        Input in{5,
                 {
                     {1, 1, 2}, {1, 2, 3}, {1, 3, 4}, {1, 4, 5},
                     {3, 1, 5},  // whole chain -> 1
                     {2, 3, 4},  // cut the middle, no cycle to save it
                     {3, 1, 5},  // -> 0
                     {3, 1, 3},  // left -> 1
                     {3, 4, 5},  // right -> 1
                     {3, 2, 4},  // straddles the cut -> 0
                 }};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 0, 1, 1, 0}));
        expect_agree(in);
    }

    // ---- edges enumerated --------------------------------------------------
    // query before any edge exists.
    {
        Input in{3, {{3, 1, 2}}};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{0}));
        expect_agree(in);
    }
    // self-connectivity: u == v is always 1, edge or no edge.
    {
        Input in{3, {{3, 2, 2}, {1, 1, 2}, {3, 1, 1}}};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 1}));
        expect_agree(in);
    }
    // churn on one edge: add/remove the same pair over and over, query between.
    {
        Input in{2,
                 {
                     {1, 1, 2}, {3, 1, 2},  // 1
                     {2, 1, 2}, {3, 1, 2},  // 0
                     {1, 1, 2}, {3, 1, 2},  // 1
                     {2, 1, 2}, {3, 1, 2},  // 0
                     {1, 1, 2}, {3, 1, 2},  // 1
                 }};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 0, 1, 0, 1}));
        expect_agree(in);
    }
    // a re-added edge that reconnects: cut a bridge, no replacement, then re-add.
    {
        Input in{3,
                 {
                     {1, 1, 2}, {1, 2, 3},
                     {3, 1, 3},  // 1
                     {2, 2, 3},  // split off 3
                     {3, 1, 3},  // 0
                     {1, 1, 3},  // re-bridge, different edge
                     {3, 1, 3},  // 1
                 }};
        CPF_CHECK((run_hdlt(in) == std::vector<int>{1, 0, 1}));
        expect_agree(in);
    }

    // ---- randomized differential vs the brute recompute oracle -------------
    // thousands of valid online programs. dense inserts build cycles, heavy deletes
    // cut tree edges and force the replacement search across levels; every type-3
    // answer is compared. one disagreement prints the seed that replays it.
    {
        constexpr int kSeqs = 12000;
        bool ok = cpf::differential(
            kSeqs, 0x46DC01,
            [](cpf::Rng& rng) { return gen(rng); },
            [](const Input& in) { return run_hdlt(in); },
            [](const Input& in) { return run_ref(in); });
        CPF_CHECK(ok);
        std::printf("HDLT differential: %d valid online programs vs brute oracle\n",
                    kSeqs);
    }

    return cpf::report();
}
