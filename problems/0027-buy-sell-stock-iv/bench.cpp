#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// n = 1e6. two regimes:
//   1) K = n/2 on random prices -- K clears the trade count, so solve peels the
//      peak in one O(n) pass and returns. the "easy" large-K case.
//   2) the real worst case: a zigzag with ~n/2 trades available and K set BELOW
//      that, forcing the full binary search -- ~log2(maxPrice) penalized sweeps
//      over the whole array. this is the O(n log maxPrice) the trick is named for.
// solve is a single heavy call, not a hot micro-loop, so time it directly with a
// steady clock and print milliseconds.
namespace {

void time_solve(const char* label, int K, const std::vector<int>& p) {
    // warm once -- first touch pays page faults and cold cache.
    volatile p0027::i64 sink = p0027::solve(K, p);
    (void)sink;

    auto t0 = std::chrono::steady_clock::now();
    p0027::i64 ans = p0027::solve(K, p);
    auto t1 = std::chrono::steady_clock::now();

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    cpf::keep(ans);
    std::printf("%-32s n=%d K=%d  ans=%lld  %8.2f ms\n", label,
                static_cast<int>(p.size()), K, static_cast<long long>(ans), ms);
}

}  // namespace

int main() {
    constexpr int kN = 1'000'000;
    cpf::Rng rng(2027);

    // regime 1: random prices, K = n/2. large K, returns after one sweep.
    {
        std::vector<int> p(kN);
        for (auto& x : p) x = static_cast<int>(rng.in_range(1, 1000000000));
        time_solve("random, K=n/2 (peak peel)", kN / 2, p);
    }

    // regime 2: zigzag -- ~n/2 isolated trades. K = n/4 sits well below that, so
    // K < r and the binary search runs its full ~30 sweeps. the honest worst case.
    {
        std::vector<int> p(kN);
        for (int i = 0; i < kN; ++i) p[i] = (i & 1) ? 1000000000 : 1;
        time_solve("zigzag, K=n/4 (binary search)", kN / 4, p);
    }

    // regime 3: same zigzag with random amplitudes, K = n/8 -- distinct slopes so
    // the search can't collapse early. worst of the worst.
    {
        std::vector<int> p(kN);
        for (int i = 0; i < kN; ++i)
            p[i] = (i & 1) ? static_cast<int>(rng.in_range(500000000, 1000000000))
                           : static_cast<int>(rng.in_range(1, 500000000));
        time_solve("random zigzag, K=n/8", kN / 8, p);
    }

    return 0;
}
