#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

using p0055::u64;

// 500 hard semiprimes: two primes near 1e9, so each product sits near 1e18 with
// no small factor to strip and one big-ish prime to find -- rho's worst case,
// ~n^{1/4} ~ 3e4 iterations to the smaller factor. this is the number that
// matters: the wall time to clear a full q = 500 batch of the ugliest inputs.
int main() {
    cpf::Rng rng(2025);
    constexpr int kN = 500;

    auto prime_near = [&](u64 lo, u64 hi) {
        while (true) {
            u64 x = static_cast<u64>(rng.in_range(
                        static_cast<std::int64_t>(lo),
                        static_cast<std::int64_t>(hi))) |
                    1ull;
            if (p0055::is_prime(x)) return x;
        }
    };

    std::vector<u64> nums;
    nums.reserve(kN);
    for (int i = 0; i < kN; ++i) {
        u64 p = prime_near(900000000ull, 1000000000ull);
        u64 q = prime_near(900000000ull, 1000000000ull);
        nums.push_back(p * q);
    }

    // batch wall time -- the answer to "how long does the hardest q = 500 take".
    auto t0 = std::chrono::steady_clock::now();
    u64 sink = 0;
    for (u64 n : nums) {
        auto f = p0055::factorize(n);
        sink += f.size();
        for (u64 x : f) sink += x;
    }
    auto t1 = std::chrono::steady_clock::now();
    cpf::keep(sink);

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("factor 500 hard ~1e18 semiprimes: %.2f ms total, %.1f us each\n",
                ms, ms * 1000.0 / kN);

    // steady-state throughput, cycling the same batch so generation stays out of
    // the timed loop.
    std::size_t i = 0;
    cpf::bench("factorize ~1e18 semiprime", 20000, [&] {
        u64 n = nums[i++ % kN];
        auto f = p0055::factorize(n);
        return f.empty() ? 0ull : f[0];
    });
    return 0;
}
