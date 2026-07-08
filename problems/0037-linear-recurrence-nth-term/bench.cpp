#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

using p0037::i64;
using p0037::kMod;

namespace {

// emit m terms of an order-k recurrence -- gives BM a real order-k law to find.
std::vector<i64> emit(const std::vector<i64>& rec, const std::vector<i64>& init,
                      int m) {
    std::vector<i64> a(init);
    a.resize(static_cast<std::size_t>(m));
    const int k = static_cast<int>(rec.size());
    for (int i = k; i < m; ++i) {
        __int128 acc = 0;
        for (int j = 0; j < k; ++j)
            acc += static_cast<__int128>(rec[j]) * a[i - 1 - j];
        a[static_cast<std::size_t>(i)] = static_cast<i64>(acc % kMod);
    }
    return a;
}

// one worst-ish case: m terms of a genuine order-k recurrence, then a_N at the
// far end of the index range. reports wall time for the full BM + kitamasa run.
void run(int m, int k, std::uint64_t N) {
    cpf::Rng rng(7);
    std::vector<i64> rec(k), init(k);
    for (int j = 0; j < k; ++j) rec[j] = rng.in_range(0, kMod - 1);
    if (rec[k - 1] == 0) rec[k - 1] = 1;
    for (int j = 0; j < k; ++j) init[j] = rng.in_range(0, kMod - 1);
    std::vector<i64> a = emit(rec, init, m);

    // one iter -- a single solve at k=few-thousand is already seconds; the bench
    // harness warms once then times once, which is what we want here.
    cpf::bench("BM+kitamasa m/k/N", 1, [&] { return p0037::solve(N, a); });
    std::printf("  (m=%d, order k=%d, N=%llu)\n", m, k,
                static_cast<unsigned long long>(N));
}

}  // namespace

int main() {
    constexpr std::uint64_t N = 1000000000000000000ULL;  // 1e18.

    // headline: the full m=2e4 read, a large recovered order, the hardest index.
    // BM is O(m*k); kitamasa is O(k^2 log N) and dominates as k climbs.
    run(20000, 2000, N);
    run(20000, 5000, N);
    run(20000, 10000, N);
    return 0;
}
