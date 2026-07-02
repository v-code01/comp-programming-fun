#include <array>
#include <cstdint>
#include <cstdio>

#include "solution.hpp"

// stdin: line 1 is W, line 2 is cnt_1..cnt_8. stdout: the max good total.
int main() {
    std::int64_t W;
    if (std::scanf("%lld", &W) != 1) return 1;
    std::array<std::int64_t, 9> cnt{};  // index 0 unused -- weights are 1..8.
    for (int i = 1; i <= 8; ++i)
        if (std::scanf("%lld", &cnt[i]) != 1) return 1;
    std::printf("%lld\n", static_cast<long long>(p0005::solve(W, cnt)));
    return 0;
}
