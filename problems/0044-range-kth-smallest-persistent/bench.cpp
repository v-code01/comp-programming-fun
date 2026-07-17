#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case is n = q = 1e5. build is n persistent inserts (the O(n log n) half);
// the query stream is q parallel root-to-leaf descents. time the two separately so
// the numbers don't blur into one another. cpf::keep pins each answer so the
// optimizer can't delete the query loop.
namespace {

constexpr int kN = 100000;
constexpr int kQ = 100000;

struct Query {
    int l;
    int r;
    int k;
};

// full-magnitude signed values, mostly distinct -- the deepest domain the tree
// will build, ~17 levels. that's the stressing shape for the insert path.
std::vector<std::int64_t> make_array(cpf::Rng& rng) {
    std::vector<std::int64_t> a(static_cast<std::size_t>(kN));
    for (auto& x : a) x = rng.in_range(-1000000000, 1000000000);
    return a;
}

// random valid queries: l <= r, 1 <= k <= r-l+1. spread across the array so ranges
// run from single elements to the whole thing.
std::vector<Query> make_queries(cpf::Rng& rng) {
    std::vector<Query> qs;
    qs.reserve(static_cast<std::size_t>(kQ));
    for (int i = 0; i < kQ; ++i) {
        const int l = static_cast<int>(rng.in_range(1, kN));
        const int r = static_cast<int>(rng.in_range(l, kN));
        const int k = static_cast<int>(rng.in_range(1, r - l + 1));
        qs.push_back({l, r, k});
    }
    return qs;
}

}  // namespace

int main() {
    cpf::Rng rng(2044);
    const std::vector<std::int64_t> a = make_array(rng);
    const std::vector<Query> qs = make_queries(rng);

    double build_ms = 1e18;
    double query_ms = 1e18;
    // best of three, rebuilding each round to shed scheduler noise.
    for (int round = 0; round < 3; ++round) {
        p0044::PersistentKth tree;

        auto b0 = std::chrono::steady_clock::now();
        tree.build(a);
        auto b1 = std::chrono::steady_clock::now();
        double bms = std::chrono::duration<double, std::milli>(b1 - b0).count();
        if (bms < build_ms) build_ms = bms;

        auto q0 = std::chrono::steady_clock::now();
        for (const Query& q : qs) cpf::keep(tree.query(q.l, q.r, q.k));
        auto q1 = std::chrono::steady_clock::now();
        double qms = std::chrono::duration<double, std::milli>(q1 - q0).count();
        if (qms < query_ms) query_ms = qms;
    }

    const double build_per = build_ms * 1e6 / static_cast<double>(kN);   // ns/insert
    const double query_per = query_ms * 1e6 / static_cast<double>(kQ);   // ns/query
    std::printf("n=%d  build %8.2f ms  %7.1f ns/insert\n", kN, build_ms, build_per);
    std::printf("q=%d  query %8.2f ms  %7.1f ns/query\n", kQ, query_ms, query_per);
    return 0;
}
