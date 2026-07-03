#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Prereqs = std::vector<std::pair<int, int>>;

// one random instance: a small project count, profits on both sides of zero, and
// a random pile of prerequisite edges. edges are drawn freely, so cycles show up
// on their own -- and a cycle is the interesting case: it welds a group into
// all-or-nothing, exactly what forces the min-cut to move whole components.
struct Instance {
    std::vector<std::int64_t> profits;
    Prereqs prereqs;
};

Instance random_instance(cpf::Rng& rng, int max_n, std::int64_t plo,
                         std::int64_t phi) {
    const int n = static_cast<int>(rng.in_range(1, max_n));
    Instance in;
    in.profits.resize(static_cast<std::size_t>(n));
    for (auto& p : in.profits) p = rng.in_range(plo, phi);

    // up to ~2n edges. self-loops (i->i) are trivially closed and harmless; we
    // skip them so the edges actually constrain something.
    const int m = static_cast<int>(rng.in_range(0, 2 * n));
    for (int e = 0; e < m; ++e) {
        const int i = static_cast<int>(rng.in_range(0, n - 1));
        const int j = static_cast<int>(rng.in_range(0, n - 1));
        if (i != j) in.prereqs.emplace_back(i, j);
    }
    return in;
}

}  // namespace

int main() {
    using namespace p0021;

    // ---- hand shapes, every "expected" cross-checked against the 2^n oracle so
    // the number is never just asserted from my head. ----

    // n = 1, single positive project, no deps -> take it.
    {
        std::vector<std::int64_t> p = {7};
        Prereqs pre = {};
        CPF_EQ(max_profit(p, pre), 7);
        CPF_EQ(max_profit_brute(p, pre), 7);
    }
    // n = 1, single negative project -> ship nothing, net 0.
    {
        std::vector<std::int64_t> p = {-7};
        Prereqs pre = {};
        CPF_EQ(max_profit(p, pre), 0);
        CPF_EQ(max_profit_brute(p, pre), 0);
    }

    // all-negative -> the empty closure wins, answer 0.
    {
        std::vector<std::int64_t> p = {-3, -1, -8, -2, -5};
        Prereqs pre = {{0, 1}, {2, 3}, {3, 4}};
        CPF_EQ(max_profit(p, pre), 0);
        CPF_EQ(max_profit_brute(p, pre), 0);
    }

    // all-positive, no deps -> take everything, answer is the sum.
    {
        std::vector<std::int64_t> p = {3, 1, 8, 2, 5};
        Prereqs pre = {};
        CPF_EQ(max_profit(p, pre), 19);
        CPF_EQ(max_profit_brute(p, pre), 19);
    }

    // a chain 0->1->2->3->4 of mixed signs. selecting 0 drags the whole tail in.
    // profits 10, -4, 3, -1, -2. taking {0..4} nets 10-4+3-1-2 = 6. can you stop
    // earlier? 0 alone needs 1,2,3,4. 2 alone needs 3,4 -> 3-1-2 = 0. so 6 is it,
    // and it beats the empty 0. the oracle confirms.
    {
        std::vector<std::int64_t> p = {10, -4, 3, -1, -2};
        Prereqs pre = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
        CPF_EQ(max_profit_brute(p, pre), 6);
        CPF_EQ(max_profit(p, pre), 6);
    }

    // chain where the tail is too expensive -> ship nothing. 0->1->2, profits
    // 5, -4, -6: taking 0 forces 1 and 2 -> 5-4-6 = -5 < 0. answer 0.
    {
        std::vector<std::int64_t> p = {5, -4, -6};
        Prereqs pre = {{0, 1}, {1, 2}};
        CPF_EQ(max_profit_brute(p, pre), 0);
        CPF_EQ(max_profit(p, pre), 0);
    }

    // a prerequisite CYCLE: 0->1->2->0 welds {0,1,2} into all-or-nothing. profits
    // 10, -3, -2 sum to 5 > 0, so the whole group ships -> 5.
    {
        std::vector<std::int64_t> p = {10, -3, -2};
        Prereqs pre = {{0, 1}, {1, 2}, {2, 0}};
        CPF_EQ(max_profit_brute(p, pre), 5);
        CPF_EQ(max_profit(p, pre), 5);
    }
    // same cycle, now the group is a net loss (2, -3, -2 = -3) -> ship nothing.
    {
        std::vector<std::int64_t> p = {2, -3, -2};
        Prereqs pre = {{0, 1}, {1, 2}, {2, 0}};
        CPF_EQ(max_profit_brute(p, pre), 0);
        CPF_EQ(max_profit(p, pre), 0);
    }
    // two cycles, one profitable and independent of the other. {0,1} cycle nets
    // 9-2 = 7; {2,3} cycle nets 1-5 = -4 (dropped). answer 7.
    {
        std::vector<std::int64_t> p = {9, -2, 1, -5};
        Prereqs pre = {{0, 1}, {1, 0}, {2, 3}, {3, 2}};
        CPF_EQ(max_profit_brute(p, pre), 7);
        CPF_EQ(max_profit(p, pre), 7);
    }

    // a shared prerequisite paid once: 0 and 1 are both profitable and both need
    // 2 (a cost). 6, 4, -3: take {0,1,2} -> 6+4-3 = 7. the single -3 is amortized
    // across both, and INF edges make the cut carry it exactly once.
    {
        std::vector<std::int64_t> p = {6, 4, -3};
        Prereqs pre = {{0, 2}, {1, 2}};
        CPF_EQ(max_profit_brute(p, pre), 7);
        CPF_EQ(max_profit(p, pre), 7);
    }

    // profit exactly zero on a node -- free to include, never moves the answer.
    {
        std::vector<std::int64_t> p = {0, 5, -1};
        Prereqs pre = {{1, 2}};
        CPF_EQ(max_profit_brute(p, pre), 4);
        CPF_EQ(max_profit(p, pre), 4);
    }

    // big magnitudes: the INF-vs-P headroom. P = 5 * 1e9 sits far under 1e18.
    {
        std::vector<std::int64_t> p = {1000000000LL, 1000000000LL, 1000000000LL,
                                       1000000000LL, 1000000000LL};
        Prereqs pre = {{0, 1}, {1, 2}, {2, 3}, {3, 4}};
        CPF_EQ(max_profit(p, pre), 5000000000LL);
        CPF_EQ(max_profit_brute(p, pre), 5000000000LL);
    }
    // big mixed: one huge positive dragging one huge negative -> net still positive.
    {
        std::vector<std::int64_t> p = {1000000000LL, -999999999LL};
        Prereqs pre = {{0, 1}};
        CPF_EQ(max_profit(p, pre), 1);
        CPF_EQ(max_profit_brute(p, pre), 1);
    }

    // ---- randomized differential vs the 2^n oracle. thousands of small
    // instances, profits on both sides of zero, random edges -> cycles for free.
    // small n keeps the oracle's 2^n cheap. ----
    {
        const int kCases = 20000;
        bool ok = cpf::differential(
            kCases, 0x2121u,
            [](cpf::Rng& rng) { return random_instance(rng, 12, -20, 20); },
            [](const Instance& in) { return max_profit(in.profits, in.prereqs); },
            [](const Instance& in) {
                return max_profit_brute(in.profits, in.prereqs);
            });
        CPF_CHECK(ok);
        std::printf("randomized differential vs 2^n oracle: %d cases, n<=12\n",
                    kCases);
    }

    // ---- a second sweep skewed dense and tiny (n<=8, more edges than nodes) so
    // prerequisite cycles and forced groups dominate. count how many land on the
    // empty-set answer to prove the 0 floor is exercised, not just claimed. ----
    {
        const int kCases = 20000;
        int zero_answers = 0;
        int diffs = 0;
        cpf::Rng rng(0x515Eu);
        for (int c = 0; c < kCases; ++c) {
            // denser edges: up to 3n, tiny n, so groups weld constantly.
            const int n = static_cast<int>(rng.in_range(2, 8));
            Instance in;
            in.profits.resize(static_cast<std::size_t>(n));
            for (auto& p : in.profits) p = rng.in_range(-15, 15);
            const int m = static_cast<int>(rng.in_range(0, 3 * n));
            for (int e = 0; e < m; ++e) {
                const int i = static_cast<int>(rng.in_range(0, n - 1));
                const int j = static_cast<int>(rng.in_range(0, n - 1));
                if (i != j) in.prereqs.emplace_back(i, j);
            }
            const std::int64_t got = max_profit(in.profits, in.prereqs);
            const std::int64_t want = max_profit_brute(in.profits, in.prereqs);
            if (want == 0) ++zero_answers;
            if (got != want) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("dense/tiny differential n<=8: %d cases, %d hit the 0 floor, "
                    "%d diffs\n",
                    kCases, zero_answers, diffs);
    }

    // ---- large-magnitude differential: profits at the 1e9 scale so the int64
    // capacities and the INF>P margin get stress at real bounds, still diffed
    // against the oracle on small n. ----
    {
        const int kCases = 5000;
        bool ok = cpf::differential(
            kCases, 0x9E9Eu,
            [](cpf::Rng& rng) {
                return random_instance(rng, 10, -1000000000LL, 1000000000LL);
            },
            [](const Instance& in) { return max_profit(in.profits, in.prereqs); },
            [](const Instance& in) {
                return max_profit_brute(in.profits, in.prereqs);
            });
        CPF_CHECK(ok);
        std::printf("large-magnitude differential: %d cases, |p|<=1e9, n<=10\n",
                    kCases);
    }

    return cpf::report();
}
