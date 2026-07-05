#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 "n m"; then n lines of m chars, '.' free and '#' blocked.
// stdout: the count of hamiltonian cycles over the free cells, mod 1e9+7.
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n = 0, m = 0;
    if (!(std::cin >> n >> m)) return 1;
    if (n <= 0 || m <= 0) {
        std::printf("0\n");
        return 0;
    }

    std::vector<std::string> grid(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!(std::cin >> grid[static_cast<std::size_t>(i)])) return 1;
    }

    std::printf("%lld\n",
                static_cast<long long>(p0034::count_hamiltonian_cycles(grid)));
    return 0;
}
