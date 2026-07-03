#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

struct Tree {
    int n = 0;
    std::vector<int> color;
    std::vector<std::pair<int, int>> edges;
};

std::vector<int> colors(int n, std::int64_t alpha, cpf::Rng& rng) {
    std::vector<int> c(static_cast<std::size_t>(n));
    for (auto& x : c) x = static_cast<int>(rng.in_range(1, alpha));
    return c;
}

// path: depth == n. the point isn't runtime -- every edge is heavy, so light
// work is nil -- it's proving the iterative solve survives a 1e5-deep chain where
// a recursive dfs would blow the stack.
Tree make_path(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t{n, colors(n, alpha, rng), {}};
    for (int v = 2; v <= n; ++v) t.edges.emplace_back(v - 1, v);
    return t;
}

// complete binary tree: parent(v) = v/2. every vertex sits under ~log n light
// edges, so the small-to-large charge is fully loaded -- this is the Theta(n log
// n) worst case the lower bound describes.
Tree make_binary(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t{n, colors(n, alpha, rng), {}};
    for (int v = 2; v <= n; ++v) t.edges.emplace_back(v / 2, v);
    return t;
}

// star: one O(n) merge at the root, everything else a leaf. the cheap extreme --
// a single fat re-insert, no depth.
Tree make_star(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t{n, colors(n, alpha, rng), {}};
    for (int v = 2; v <= n; ++v) t.edges.emplace_back(1, v);
    return t;
}

// random parent in [1, v-1]: mixed heavy/light edges, the shape real inputs take.
Tree make_random(int n, std::int64_t alpha, cpf::Rng& rng) {
    Tree t{n, colors(n, alpha, rng), {}};
    for (int v = 2; v <= n; ++v)
        t.edges.emplace_back(static_cast<int>(rng.in_range(1, v - 1)), v);
    return t;
}

void time_shape(const char* name, const Tree& t) {
    // one measured pass -- the solve allocates fresh each call, so re-timing in a
    // tight loop would benchmark the allocator, not the algorithm. warm once.
    volatile std::int64_t sink = 0;
    {
        auto warm = p0013::solve(t.n, t.color, t.edges);
        sink += warm.empty() ? 0 : warm[0];
    }
    auto t0 = std::chrono::steady_clock::now();
    auto out = p0013::solve(t.n, t.color, t.edges);
    auto t1 = std::chrono::steady_clock::now();
    sink += out.empty() ? 0 : out[0];
    cpf::keep(sink);
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("  %-18s n=%d  %10.2f ms  %8.1f Mvtx/s\n", name, t.n, ms,
                static_cast<double>(t.n) / ms / 1e3);
}

}  // namespace

int main() {
    cpf::Rng rng(600);
    constexpr int kN = 100000;
    constexpr std::int64_t kAlpha = 1000;  // enough colors to keep the max moving.

    std::printf("lomsat gelral -- dsu on tree, n=%d, alpha=%lld\n", kN,
                static_cast<long long>(kAlpha));
    time_shape("path (depth 1e5)", make_path(kN, kAlpha, rng));
    time_shape("complete binary", make_binary(kN, kAlpha, rng));
    time_shape("star", make_star(kN, kAlpha, rng));
    time_shape("random parent", make_random(kN, kAlpha, rng));

    // same complete-binary worst case, single color -- max count climbs to n and
    // every insert hits the equality branch, the hottest path in add().
    time_shape("binary, alpha=1", make_binary(kN, 1, rng));
    return 0;
}
