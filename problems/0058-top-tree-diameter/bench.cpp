#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the top tree is amortized, so a benchmark has to attack the amortization, not run
// lazy inputs. n = 1e5 vertices, then ~2e5 mixed ops on top of the build.
//
//   build -- link a random tree biased DEEP: with 1-in-3 odds the new vertex hangs
//   off its predecessor, stretching long chains so expose paths are long. deep trees
//   are where expose and the rake splices actually pay their log n.
//
//   mix -- diameter queries on random vertices force full folds to the root cluster,
//   while cut+relink of a live tree edge churns the compress/rake decomposition
//   without changing the tree's shape, so every later op stays legal.
//
// timed as one destructive run (ops mutate the forest, nothing to repeat in place);
// rebuild and take the best of a few to shed scheduler noise. cpf::keep pins the
// diameters so the read path can't be elided.
namespace {

constexpr int kN = 100000;
constexpr int kMix = 200000;
constexpr std::int64_t kW = 1000000000LL;

struct Op {
    int type;  // 1 link, 2 cut, 3 diameter
    int u;
    int v;
    std::int64_t w;
};

// a random deep-ish tree, returned as its link ops plus the edge list.
std::vector<std::pair<int, int>> build_tree(cpf::Rng& rng, std::vector<Op>& ops) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(kN - 1));
    for (int i = 2; i <= kN; ++i) {
        int p = (rng.in_range(0, 2) == 0)
                    ? i - 1  // 1-in-3: extend the chain, keep it deep.
                    : static_cast<int>(rng.in_range(1, i - 1));
        std::int64_t w = rng.in_range(1, kW);
        ops.push_back({1, i, p, w});
        edges.emplace_back(i, p);
    }
    return edges;
}

// the mixed phase: ~60% diameter queries, ~40% cut+relink of a random tree edge.
std::vector<Op> mix_program(cpf::Rng& rng,
                            const std::vector<std::pair<int, int>>& edges) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kMix + kMix / 2));
    int emitted = 0;
    while (emitted < kMix) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        if (roll < 60) {
            ops.push_back({3, static_cast<int>(rng.in_range(1, kN)), 0, 0});
            ++emitted;
        } else {
            // cut then relink the same edge -- legal (split, then rejoin), shape
            // unchanged so later ops stay connected. counts as two ops.
            auto e = edges[static_cast<std::size_t>(
                rng.in_range(0, static_cast<std::int64_t>(edges.size()) - 1))];
            std::int64_t w = rng.in_range(1, kW);
            ops.push_back({2, e.first, e.second, 0});
            ops.push_back({1, e.first, e.second, w});
            emitted += 2;
        }
    }
    return ops;
}

double time_once(const std::vector<Op>& build_ops, const std::vector<Op>& mix_ops) {
    p0058::DynamicDiameter dd(kN);
    auto run = [&](const std::vector<Op>& ops) {
        for (const Op& o : ops) {
            if (o.type == 1)
                dd.link(o.u, o.v, o.w);
            else if (o.type == 2)
                dd.cut(o.u, o.v);
            else
                cpf::keep(dd.diameter(o.u));
        }
    };
    auto t0 = std::chrono::steady_clock::now();
    run(build_ops);
    run(mix_ops);
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

}  // namespace

int main() {
    cpf::Rng rng(2026);

    std::vector<Op> build_ops;
    std::vector<std::pair<int, int>> edges = build_tree(rng, build_ops);
    std::vector<Op> mix_ops = mix_program(rng, edges);

    std::size_t total = build_ops.size() + mix_ops.size();

    double best = 1e18;
    for (int r = 0; r < 3; ++r) {  // best of three, rebuild each time.
        double ms = time_once(build_ops, mix_ops);
        if (ms < best) best = ms;
    }
    double per = best * 1e6 / static_cast<double>(total);  // ns/op

    std::printf("top-tree diameter  n=%d  build=%zu ops  mix=%zu ops  total=%zu\n",
                kN, build_ops.size(), mix_ops.size(), total);
    std::printf("%-22s  %9.2f ms  %8.1f ns/op\n", "adversarial mix", best, per);
    return 0;
}
