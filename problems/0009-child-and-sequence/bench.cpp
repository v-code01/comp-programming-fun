#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// beats is amortized, so a benchmark has to attack the amortization, not just
// run random ops. two scenarios at n = m = 1e5:
//
//   adversarial mix -- the honest worst case. point-sets keep slamming cells back
//   to 1e9, each one refunding O(log V) of mod potential, while range-mods over
//   random windows keep spending it. this is the workload the potential argument
//   is actually paying for; sums are interleaved to keep the read path live.
//
//   mod storm -- full-range mods by a shrinking x. every pass forces real
//   reductions across the array before values settle, stressing the descent.
//
// each scenario is timed as a single destructive run (ops mutate the tree, so
// there's nothing to repeat in place); we rebuild and take the best of a few
// runs to shed scheduler noise. cpf::keep pins the sums so nothing is elided.
namespace {

constexpr int kN = 100000;
constexpr int kM = 100000;
constexpr std::int64_t kV = 1000000000LL;

struct Op {
    int type;
    int l;
    int r;
    std::int64_t x;
};

std::vector<std::int64_t> make_array(cpf::Rng& rng) {
    std::vector<std::int64_t> a(static_cast<std::size_t>(kN));
    for (auto& v : a) v = rng.in_range(1, kV);
    return a;
}

// ~30% sets (repopulate high), ~45% mods, ~25% sums.
std::vector<Op> mix_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kM));
    for (int q = 0; q < kM; ++q) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        int i = static_cast<int>(rng.in_range(0, kN - 1));
        int j = static_cast<int>(rng.in_range(0, kN - 1));
        int lo = i < j ? i : j;
        int hi = i < j ? j : i;
        if (roll < 30)
            ops.push_back({3, lo, 0, rng.in_range(1, kV)});  // set high
        else if (roll < 75)
            ops.push_back({2, lo, hi, rng.in_range(2, kV)});  // mod
        else
            ops.push_back({1, lo, hi, 0});  // sum
    }
    return ops;
}

// full-range mods by a descending x, forcing genuine reductions each pass, with
// periodic sums so the aggregate stays observable.
std::vector<Op> storm_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kM));
    std::int64_t x = kV;
    for (int q = 0; q < kM; ++q) {
        if (q % 8 == 0) {
            ops.push_back({1, 0, kN - 1, 0});  // whole-array sum
        } else {
            x = x > 3 ? x - x / 3 : 2;  // shrink x geometrically, floor at 2
            ops.push_back({2, 0, kN - 1, x});
            if (q % 500 == 0) x = kV;  // reset so the storm never fully dies
            (void)rng;
        }
    }
    return ops;
}

double time_once(const std::vector<std::int64_t>& init,
                 const std::vector<Op>& ops) {
    p0009::Beats tree;
    tree.build(init);
    auto t0 = std::chrono::steady_clock::now();
    for (const Op& o : ops) {
        if (o.type == 1)
            cpf::keep(tree.range_sum(o.l, o.r));
        else if (o.type == 2)
            tree.range_mod(o.l, o.r, o.x);
        else
            tree.point_set(o.l, o.x);
    }
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void report(const char* name, const std::vector<std::int64_t>& init,
            const std::vector<Op>& ops) {
    double best = 1e18;
    for (int r = 0; r < 3; ++r) {  // best of three, rebuild each time.
        double ms = time_once(init, ops);
        if (ms < best) best = ms;
    }
    double per = best * 1e6 / static_cast<double>(ops.size());  // ns/op
    std::printf("%-22s  n=m=%d  %9.2f ms  %8.1f ns/op\n", name, kM, best, per);
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    std::vector<std::int64_t> init = make_array(rng);

    std::vector<Op> mix = mix_program(rng);
    std::vector<Op> storm = storm_program(rng);

    report("adversarial mix", init, mix);
    report("mod storm", init, storm);
    return 0;
}
