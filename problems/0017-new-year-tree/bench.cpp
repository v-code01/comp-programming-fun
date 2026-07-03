#include <chrono>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case is n = m = 4e5. two shapes, because they stress different halves:
//
//   path -- vertex i hangs off i-1, so the tree is one line of depth 4e5. this is
//   the shape that would overflow a recursive Euler tour; the iterative tour has
//   to swallow it flat. subtree windows here range from the whole array down to a
//   single cell, so the segment tree walks its full height constantly.
//
//   random tree -- vertex i hangs off a uniform earlier vertex, a bushy shallow
//   shape. windows are scattered and small, exercising the repaint/query mix on
//   arbitrary sub-ranges rather than nested prefixes.
//
// each op stream is ~2/3 repaints, 1/3 queries, over a 60-color alphabet so masks
// stay full. timed as one destructive run (ops mutate the tree); best of three,
// rebuilding each time, to shed scheduler noise. cpf::keep pins the query results
// so nothing is elided.
namespace {

constexpr int kN = 400000;
constexpr int kM = 400000;

struct Op {
    int type;  // 1 = repaint, 2 = distinct
    int v;
    int c;
};

std::vector<int> make_colors(cpf::Rng& rng) {
    std::vector<int> color(static_cast<std::size_t>(kN) + 1, 0);
    for (int i = 1; i <= kN; ++i)
        color[static_cast<std::size_t>(i)] = static_cast<int>(rng.in_range(1, 60));
    return color;
}

std::vector<std::pair<int, int>> path_tree() {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 2; i <= kN; ++i) edges.emplace_back(i - 1, i);
    return edges;
}

std::vector<std::pair<int, int>> random_tree(cpf::Rng& rng) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 2; i <= kN; ++i)
        edges.emplace_back(static_cast<int>(rng.in_range(1, i - 1)), i);
    return edges;
}

std::vector<Op> program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kM));
    for (int q = 0; q < kM; ++q) {
        int v = static_cast<int>(rng.in_range(1, kN));
        if (rng.in_range(0, 2) == 0)
            ops.push_back({2, v, 0});
        else
            ops.push_back({1, v, static_cast<int>(rng.in_range(1, 60))});
    }
    return ops;
}

double time_once(int n, const std::vector<int>& color,
                 const std::vector<std::pair<int, int>>& edges,
                 const std::vector<Op>& ops) {
    p0017::Tree tree;
    auto t0 = std::chrono::steady_clock::now();
    tree.build(n, color, edges);  // time the iterative tour + segtree build too.
    for (const Op& o : ops) {
        if (o.type == 1)
            tree.repaint(o.v, o.c);
        else
            cpf::keep(tree.distinct(o.v));
    }
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void report(const char* name, int n, const std::vector<int>& color,
            const std::vector<std::pair<int, int>>& edges,
            const std::vector<Op>& ops) {
    double best = 1e18;
    for (int r = 0; r < 3; ++r) {
        double ms = time_once(n, color, edges, ops);
        if (ms < best) best = ms;
    }
    double per = best * 1e6 / static_cast<double>(ops.size());  // ns/op
    std::printf("%-16s  n=m=%d  %9.2f ms  %8.1f ns/op\n", name, kM, best, per);
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    std::vector<int> color = make_colors(rng);
    std::vector<Op> ops = program(rng);

    std::vector<std::pair<int, int>> path = path_tree();
    std::vector<std::pair<int, int>> rand_t = random_tree(rng);

    report("path (deep)", kN, color, path, ops);
    report("random (bushy)", kN, color, rand_t, ops);
    return 0;
}
