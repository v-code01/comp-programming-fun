#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/rng.hpp"
#include "solution.hpp"

// the method is amortized over the whole timeline, so the benchmark attacks the
// timeline: n = q = 2e5, and the ops churn a fixed pool of edges on and off. every
// toggle shatters an edge's life into another short interval, so the pool sprays
// tens of thousands of intervals across the segtree -- the O(q log q) placement
// load this structure exists to eat. queries between random vertices land at every
// depth of the DFS. timed as one offline solve; best of three to shed noise.
namespace {

constexpr int kN = 200000;
constexpr int kQ = 200000;
constexpr int kPool = 200000;  // candidate edges we flip on and off.

using p0039::Op;

// build a churny, always-valid program. a pool edge that's off gets added, an edge
// that's on gets removed -- so type-2 never touches a dead edge. ~50% toggles,
// ~50% queries.
std::vector<Op> build_program(cpf::Rng& rng, std::vector<std::pair<int, int>>& pool) {
    pool.resize(static_cast<std::size_t>(kPool));
    for (auto& e : pool) {
        int a = static_cast<int>(rng.in_range(1, kN));
        int b = static_cast<int>(rng.in_range(1, kN));
        if (a == b) b = a % kN + 1;  // no self-loops in the pool.
        e = {a, b};
    }
    std::vector<char> on(static_cast<std::size_t>(kPool), 0);

    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kQ));
    for (int i = 0; i < kQ; ++i) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 50) {
            int idx = static_cast<int>(rng.in_range(0, kPool - 1));
            auto e = pool[static_cast<std::size_t>(idx)];
            if (on[static_cast<std::size_t>(idx)]) {
                ops.push_back({2, e.first, e.second});  // live -> remove.
                on[static_cast<std::size_t>(idx)] = 0;
            } else {
                ops.push_back({1, e.first, e.second});  // dead -> add.
                on[static_cast<std::size_t>(idx)] = 1;
            }
        } else {
            int u = static_cast<int>(rng.in_range(1, kN));
            int v = static_cast<int>(rng.in_range(1, kN));
            ops.push_back({3, u, v});
        }
    }
    return ops;
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    std::vector<std::pair<int, int>> pool;
    std::vector<Op> ops = build_program(rng, pool);

    int adds = 0, removes = 0, queries = 0;
    for (const Op& o : ops) {
        if (o.type == 1) ++adds;
        else if (o.type == 2) ++removes;
        else ++queries;
    }

    double best = 1e18;
    std::size_t out_size = 0;
    for (int r = 0; r < 3; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        std::vector<int> ans = p0039::solve(kN, ops);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best) best = ms;
        out_size = ans.size();
    }
    double per = best * 1e6 / static_cast<double>(kQ);  // ns/op

    std::printf("offline dynamic connectivity  n=%d  q=%d\n", kN, kQ);
    std::printf("  ops: %d adds  %d removes  %d queries\n", adds, removes, queries);
    std::printf("%-24s  %9.2f ms  %8.1f ns/op  (%zu answers)\n",
                "adversarial churn", best, per, out_size);
    return 0;
}
