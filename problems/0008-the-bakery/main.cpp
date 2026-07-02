#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: line 1 is "n k", line 2 is n cake types. stdout: the max total value.
int main() {
    int n, k;
    if (std::scanf("%d %d", &n, &k) != 2) return 1;
    std::vector<int> a(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (std::scanf("%d", &a[static_cast<std::size_t>(i)]) != 1) return 1;
    }
    std::printf("%d\n", p0008::solve(a, k));
    return 0;
}
