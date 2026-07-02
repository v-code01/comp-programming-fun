#include <cstddef>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: line 1 is "n k", line 2 is n prices. stdout: every achievable total for
// exactly k picks, ascending, space-separated on one line.
int main() {
    int n, k;
    if (std::scanf("%d %d", &n, &k) != 2) return 1;
    std::vector<int> a(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        if (std::scanf("%d", &a[static_cast<std::size_t>(i)]) != 1) return 1;

    std::vector<int> out = p0010::solve(k, a);
    for (std::size_t i = 0; i < out.size(); ++i)
        std::printf("%d%c", out[i], i + 1 == out.size() ? '\n' : ' ');
    return 0;
}
