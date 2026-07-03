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
//   path -- the deepest HLD chains, and the shape where an accidental recursion
//   in the build would blow the stack. it builds and runs, so it doubles as the
//   stack-safety proof: a 1e5-vertex path survives.
//
//   random tree -- typical shape, many short chains, the common query path.
//
// build is timed on its own (the O(n) half), then a 1e5-op program of interleaved
// paints and queries is timed as one run. cpf::keep pins query answers so the
// optimizer can't delete the query loop. paints mutate the tree, so the program
// runs once per shape -- nothing to repeat in place.
namespace {

constexpr int kN = 100000;
constexpr int kM = 100000;

struct Op {
    int type;
    int a;
    int b;
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

// ~40% paints (both colors), ~60% path queries -- every op walks the same
// O(log^2 n) chain machinery, queries the heavier of the two.
std::vector<Op> mix_program(cpf::Rng& rng) {
    const int ne = kN - 1;
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kM));
    for (int q = 0; q < kM; ++q) {
        const int roll = static_cast<int>(rng.in_range(0, 4));
        if (roll < 2) {  // paint.
            const int type = rng.in_range(0, 1) == 0 ? 1 : 2;
            const int e = static_cast<int>(rng.in_range(0, ne - 1));
            ops.push_back({type, e, 0});
        } else {  // query.
            const int a = static_cast<int>(rng.in_range(0, kN - 1));
            const int b = static_cast<int>(rng.in_range(0, kN - 1));
            ops.push_back({3, a, b});
        }
    }
    return ops;
}

void report(const char* name, const std::vector<std::pair<int, int>>& edges,
            const std::vector<Op>& ops) {
    double build_ms = 1e18;
    double ops_ms = 1e18;
    // best of three, rebuilding each time to shed scheduler noise.
    for (int r = 0; r < 3; ++r) {
        p0018::BeardHLD tree;

        auto b0 = std::chrono::steady_clock::now();
        tree.build(kN, edges);
        auto b1 = std::chrono::steady_clock::now();
        double bms = std::chrono::duration<double, std::milli>(b1 - b0).count();
        if (bms < build_ms) build_ms = bms;

        auto o0 = std::chrono::steady_clock::now();
        for (const Op& o : ops) {
            if (o.type == 1)
                tree.paint_black(o.a);
            else if (o.type == 2)
                tree.paint_white(o.a);
            else
                cpf::keep(tree.query(o.a, o.b));
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
