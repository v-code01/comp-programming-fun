#include <cstdint>
#include <cstdio>

#include "solution.hpp"

// stdin: four lines -- counts of "((", "()", ")(", "))". stdout: 1 or 0.
int main() {
    std::int64_t c1, c2, c3, c4;
    if (std::scanf("%lld %lld %lld %lld", &c1, &c2, &c3, &c4) != 4) return 1;
    std::printf("%d\n", p0001::regular(c1, c2, c3, c4));
    return 0;
}
