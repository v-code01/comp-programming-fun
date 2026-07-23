#include <cstdio>

#include "solution.hpp"

// stdin: q, then q integers n (1 <= n <= 1e18). stdout: one line per n -- its
// prime factors, non-decreasing, space-separated. n == 1 prints an empty line.
int main() {
    int q;
    if (std::scanf("%d", &q) != 1) return 1;
    for (int i = 0; i < q; ++i) {
        unsigned long long n;
        if (std::scanf("%llu", &n) != 1) return 1;
        auto f = p0055::factorize(n);
        for (std::size_t j = 0; j < f.size(); ++j) {
            if (j) std::putchar(' ');
            std::printf("%llu", static_cast<unsigned long long>(f[j]));
        }
        std::putchar('\n');
    }
    return 0;
}
