#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// beats is amortized, so a benchmark has to attack the amortization, not just
// run random ops. three scenarios at n = q = 1e5:
//
//   chmin storm -- the honest adversary. adds keep lifting a fresh maximum band
//   while range-chmins keep shaving it, so the potential the beat argument pays
//   for is spent and refunded on every pass. this is the workload that actually
//   exercises the log^2 term; sums and maxes are interleaved to keep the read
//   path live and nothing elided.
//
//   descending cap -- full-range chmins by a shrinking x. each pass forces the
//   top band down before it settles, driving genuine se-level recursion instead
//   of cheap O(1) beats.
//
//   mixed -- the balanced case: a quarter each of chmin, add, max, sum over
//   random windows, values kept in a middling band so ties are frequent.
//
// each scenario is a single destructive run (ops mutate the tree), rebuilt and
// taken best-of-three to shed scheduler noise. cpf::keep pins the answers so the
// optimizer can't delete the read path.
namespace {

constexpr int kN = 100000;
constexpr int kQ = 100000;
constexpr std::int64_t kV = 1000000000LL;

struct Op {
    int type;
    int l;
    int r;
    std::int64_t x;
};

std::vector<std::int64_t> make_array(cpf::Rng& rng) {
    std::vector<std::int64_t> a(static_cast<std::size_t>(kN));
    for (auto& v : a) v = rng.in_range(-kV, kV);
    return a;
}

// ~30% adds (lift a band), ~45% chmins (shave it), ~15% max, ~10% sum. the adds
// and chmins over random windows are what keep potential churning.
std::vector<Op> storm_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kQ));
    for (int q = 0; q < kQ; ++q) {
        int roll = static_cast<int>(rng.in_range(0, 99));
        int i = static_cast<int>(rng.in_range(0, kN - 1));
        int j = static_cast<int>(rng.in_range(0, kN - 1));
        int lo = i < j ? i : j;
        int hi = i < j ? j : i;
        if (roll < 30)
            ops.push_back({2, lo, hi, rng.in_range(1, 1000)});  // add a band up
        else if (roll < 75)
            ops.push_back({1, lo, hi, rng.in_range(-kV, kV)});  // chmin it down
        else if (roll < 90)
            ops.push_back({3, lo, hi, 0});  // max
        else
            ops.push_back({4, lo, hi, 0});  // sum
    }
    return ops;
}

// full-range chmins by a descending x, forcing real reductions each pass, with
// periodic whole-array sums so the aggregate stays observable.
std::vector<Op> descending_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kQ));
    std::int64_t x = kV;
    for (int q = 0; q < kQ; ++q) {
        if (q % 8 == 0) {
            ops.push_back({4, 0, kN - 1, 0});  // whole-array sum
        } else {
            x = x > -kV ? x - kV / 200 : kV;  // slide the cap down, then reset
            ops.push_back({1, 0, kN - 1, x});
            (void)rng;
        }
    }
    return ops;
}

// balanced quarters over random windows, values in a middling band.
std::vector<Op> mixed_program(cpf::Rng& rng) {
    std::vector<Op> ops;
    ops.reserve(static_cast<std::size_t>(kQ));
    for (int q = 0; q < kQ; ++q) {
        int t = static_cast<int>(rng.in_range(1, 4));
        int i = static_cast<int>(rng.in_range(0, kN - 1));
        int j = static_cast<int>(rng.in_range(0, kN - 1));
        int lo = i < j ? i : j;
        int hi = i < j ? j : i;
        std::int64_t x = 0;
        if (t == 1)
            x = rng.in_range(-100000, 100000);
        else if (t == 2)
            x = rng.in_range(-1000, 1000);
        ops.push_back({t, lo, hi, x});
    }
    return ops;
}

double time_once(const std::vector<std::int64_t>& init,
                 const std::vector<Op>& ops) {
    p0045::Beats tree;
    tree.build(init);
    auto t0 = std::chrono::steady_clock::now();
    for (const Op& o : ops) {
        if (o.type == 1)
            tree.chmin(o.l, o.r, o.x);
        else if (o.type == 2)
            tree.add(o.l, o.r, o.x);
        else if (o.type == 3)
            cpf::keep(tree.query_max(o.l, o.r));
        else
            cpf::keep(tree.query_sum(o.l, o.r));
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
    std::printf("%-18s  n=q=%d  %9.2f ms  %8.1f ns/op\n", name, kQ, best, per);
}

}  // namespace

int main() {
    cpf::Rng rng(2026);
    std::vector<std::int64_t> init = make_array(rng);

    std::vector<Op> storm = storm_program(rng);
    std::vector<Op> descending = descending_program(rng);
    std::vector<Op> mixed = mixed_program(rng);

    report("chmin storm", init, storm);
    report("descending cap", init, descending);
    report("mixed", init, mixed);
    return 0;
}
