#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case from the constraints: n = 1e5, and queries whose k sums to 1e5.
// spread it as 1000 queries of k = 100 each so both the per-query sort and the
// virtual-tree build get exercised, not one giant query. the tree is a random
// parent tree -- realistic depth, real lca work. pre-generate everything so the
// timed number is build + all queries, not the rng.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 100000;
    constexpr int kQueries = 1000;
    constexpr int kPerQuery = 100;  // sum k = 1e5.

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 2; i <= kN; ++i) {
        edges.emplace_back(static_cast<int>(rng.in_range(1, i - 1)), i);
    }

    // each query is kPerQuery distinct nodes. dupes would shrink k, so draw a
    // shuffled prefix.
    std::vector<std::vector<int>> queries(kQueries);
    std::vector<int> pool(static_cast<std::size_t>(kN));
    for (int i = 0; i < kN; ++i) pool[static_cast<std::size_t>(i)] = i + 1;
    for (auto& q : queries) {
        q.resize(kPerQuery);
        for (int j = 0; j < kPerQuery; ++j) {
            const int r = static_cast<int>(rng.in_range(j, kN - 1));
            std::swap(pool[static_cast<std::size_t>(j)], pool[static_cast<std::size_t>(r)]);
            q[static_cast<std::size_t>(j)] = pool[static_cast<std::size_t>(j)];
        }
    }

    // build once outside the timed answer loop -- build is a one-time O(n log n),
    // the queries are the hot path the constraints actually scale.
    p0016::Kingdom kingdom;
    kingdom.build(kN, edges);

    // time: answer all 1000 queries (sum k = 1e5). one iter is the full sweep.
    cpf::bench("kingdom answer sum_k=1e5", 5, [&] {
        std::int64_t acc = 0;
        for (auto& q : queries) {
            auto copy = q;  // query sorts in place; keep the pool reusable.
            acc += kingdom.query(copy);
        }
        return acc;
    });
    std::printf("(ns/op above is all %d queries; divide by 1e6 for ms)\n", kQueries);

    // and the one-time build cost on its own, for the full picture.
    cpf::bench("kingdom build n=1e5", 3, [&] {
        p0016::Kingdom k;
        k.build(kN, edges);
        return 0;
    });
    return 0;
}
