#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// HDLT is amortized, so the benchmark has to attack the amortization, not just run
// random ops. n = 1e5 vertices, q = 2e5 ops.
//
//   build -- lay down a dense-ish graph. a spanning path keeps it connected, then
//   pile on random chords so most edges sit on cycles. an edge on a cycle, when
//   deleted, forces the replacement search to actually run and push levels -- a
//   bridge would just split for free and dodge the expensive path.
//
//   churn -- delete a present edge and insert an absent one, interleaved with
//   connectivity queries. the delete of a cycle edge is the adversarial case: cut
//   the tree edge, then scan and re-level to find the replacement. every op stays
//   legal because the present-set is tracked, so deletes hit live edges and inserts
//   hit absent ones.
//
// ops are generated up front and the replay alone is timed -- the rng and set
// bookkeeping don't belong in the number. cpf::keep pins the query results so the
// read path can't be elided.
namespace {

constexpr int kN = 100000;
constexpr int kQ = 200000;

struct Op {
    int type;  // 1 insert, 2 delete, 3 query
    int u;
    int v;
};

std::uint64_t ek(int u, int v) {
    int a = u < v ? u : v, b = u < v ? v : u;
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
}

std::vector<Op> build_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kN + kQ));

    std::unordered_set<std::uint64_t> present;
    present.reserve(static_cast<std::size_t>(kN) * 6);
    std::vector<std::pair<int, int>> live;  // present edges, for delete picks.
    live.reserve(static_cast<std::size_t>(kN) * 6);

    auto add_edge = [&](int u, int v) {
        if (u == v) return false;
        std::uint64_t k = ek(u, v);
        if (present.count(k)) return false;
        present.insert(k);
        live.emplace_back(u, v);
        ops.push_back({1, u, v});
        return true;
    };

    // spanning path 1-2-...-n keeps the whole graph connected -- ~n inserts, the
    // unavoidable cost of connecting 1e5 vertices online.
    for (int i = 1; i < kN; ++i) add_edge(i, i + 1);
    // chords so most edges land on cycles, capped so the churn phase owns the bulk
    // of the op budget.
    constexpr int kChords = 30000;
    for (int i = 0; i < kChords; ++i) {
        int u = static_cast<int>(rng.in_range(1, kN));
        int v = static_cast<int>(rng.in_range(1, kN));
        add_edge(u, v);
    }

    auto drop_live = [&](std::size_t idx) {
        live[idx] = live.back();
        live.pop_back();
    };

    // churn phase: fill to kQ total ops with delete+insert pairs and queries.
    while (static_cast<int>(ops.size()) < kQ) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 40 && !live.empty()) {
            // delete a present edge -- often a cycle edge, the pricey case.
            std::size_t idx = static_cast<std::size_t>(
                rng.in_range(0, static_cast<std::int64_t>(live.size()) - 1));
            auto [u, v] = live[idx];
            present.erase(ek(u, v));
            drop_live(idx);
            ops.push_back({2, u, v});
        } else if (roll < 75) {
            // insert an absent edge back in, keeping the graph dense.
            for (int tries = 0; tries < 8; ++tries) {
                int u = static_cast<int>(rng.in_range(1, kN));
                int v = static_cast<int>(rng.in_range(1, kN));
                if (add_edge(u, v)) break;
            }
        } else {
            int u = static_cast<int>(rng.in_range(1, kN));
            int v = static_cast<int>(rng.in_range(1, kN));
            ops.push_back({3, u, v});
        }
    }
    return ops;
}

double time_once(const std::vector<Op>& ops, int n) {
    p0046::DynamicConnectivity dc;
    dc.init(n);
    auto t0 = std::chrono::steady_clock::now();
    for (const Op& o : ops) {
        if (o.type == 1)
            dc.insert(o.u, o.v);
        else if (o.type == 2)
            dc.remove(o.u, o.v);
        else
            cpf::keep(dc.connected(o.u, o.v));
    }
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    std::vector<Op> ops = build_program(rng);

    int inserts = 0, deletes = 0, queries = 0;
    for (const Op& o : ops) {
        if (o.type == 1)
            ++inserts;
        else if (o.type == 2)
            ++deletes;
        else
            ++queries;
    }

    double best = 1e18;
    for (int r = 0; r < 3; ++r) {  // best of three, fresh structure each time.
        double ms = time_once(ops, kN);
        if (ms < best) best = ms;
    }
    double per = best * 1e6 / static_cast<double>(ops.size());  // ns/op

    std::printf("online dynamic connectivity  n=%d  ops=%zu\n", kN, ops.size());
    std::printf("  inserts=%d deletes=%d queries=%d\n", inserts, deletes, queries);
    std::printf("%-24s  %9.2f ms  %8.1f ns/op\n", "adversarial churn", best, per);
    return 0;
}
