#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case is n = m = 1e5. two tree shapes stress different parts:
//
//   path -- the deepest centroid tree (height ~log n, but every level's BFS spans
//   a long chain). this is where an accidental recursion would blow the stack, so
//   it doubles as the stack-safety proof: it builds and runs, it survives.
//
//   random tree -- typical shape, shallower pieces, tests the common path.
//
// build is timed on its own (it's the O(n log n) half), then a 1e5-op program of
// interleaved paints and queries is timed as one run. cpf::keep pins query
// answers so the optimizer can't delete the query loop. ops mutate best[], so the
// program runs once per shape -- nothing to repeat in place.
namespace {

constexpr int kN = 100000;
constexpr int kM = 100000;

struct Op {
    int type;
    int v;
};

std::vector<std::pair<int, int>> path_tree() {
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 1; i < kN; ++i) e.emplace_back(i - 1, i);
    return e;
}

std::vector<std::pair<int, int>> random_tree(cpf::Rng& rng) {
    std::vector<std::pair<int, int>> e;
    e.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 1; i < kN; ++i) {
        e.emplace_back(static_cast<int>(rng.in_range(0, i - 1)), i);
    }
    return e;
}

// ~40% paints, ~60% queries -- both op kinds walk the same O(log n) ancestor list.
std::vector<Op> mix_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kM));
    for (int q = 0; q < kM; ++q) {
        const bool paint = rng.in_range(0, 4) < 2;  // 2 in 5.
        const int v = static_cast<int>(rng.in_range(0, kN - 1));
        ops.push_back({paint ? 1 : 2, v});
    }
    return ops;
}

void report(const char* name, const std::vector<std::pair<int, int>>& edges,
            const std::vector<Op>& ops) {
    double build_ms = 1e18;
    double ops_ms = 1e18;
    // best of three, rebuilding each time to shed scheduler noise.
    for (int r = 0; r < 3; ++r) {
        p0014::NearestRed tree;

        auto b0 = std::chrono::steady_clock::now();
        tree.build(kN, edges);
        auto b1 = std::chrono::steady_clock::now();
        double bms = std::chrono::duration<double, std::milli>(b1 - b0).count();
        if (bms < build_ms) build_ms = bms;

        auto o0 = std::chrono::steady_clock::now();
        for (const Op& o : ops) {
            if (o.type == 1)
                tree.paint(o.v);
            else
                cpf::keep(tree.query(o.v));
        }
        auto o1 = std::chrono::steady_clock::now();
        double oms = std::chrono::duration<double, std::milli>(o1 - o0).count();
        if (oms < ops_ms) ops_ms = oms;
    }
    const double per = ops_ms * 1e6 / static_cast<double>(ops.size());  // ns/op
    std::printf("%-14s  n=%d  build %8.2f ms   m=%d ops %8.2f ms  %7.1f ns/op\n",
                name, kN, build_ms, kM, ops_ms, per);
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    const std::vector<Op> ops = mix_program(rng);

    report("path", path_tree(), ops);
    report("random tree", random_tree(rng), ops);
    return 0;
}
