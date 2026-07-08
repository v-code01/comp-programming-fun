#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 = "m N",  line 2 = a_0 … a_{m-1}   (each already in [0, p)).
// stdout: a_N mod p.
//
// N reaches 1e18 -- under the int64 ceiling (~9.2e18), so read it signed and hand
// it to the solver as unsigned. m up to 2e4, scanf is nowhere near the bottleneck.
int main() {
    long long m = 0, N = 0;
    if (std::scanf("%lld %lld", &m, &N) != 2) return 1;

    std::vector<std::int64_t> a(static_cast<std::size_t>(m));
    for (long long i = 0; i < m; ++i) {
        long long x = 0;
        if (std::scanf("%lld", &x) != 1) return 1;
        a[static_cast<std::size_t>(i)] = x % p0037::kMod;
    }

    std::printf("%lld\n",
                static_cast<long long>(p0037::solve(static_cast<std::uint64_t>(N), a)));
    return 0;
}
