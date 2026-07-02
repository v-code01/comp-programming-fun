#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "solution.hpp"

// stdin: "n q", then q lines of "l r". stdout: max painted sections after
// firing the best two painters.
int main() {
    int n, q;
    if (std::scanf("%d %d", &n, &q) != 2) return 1;
    std::vector<std::pair<int, int>> ranges(static_cast<std::size_t>(q));
    for (int i = 0; i < q; ++i) {
        if (std::scanf("%d %d", &ranges[static_cast<std::size_t>(i)].first,
                       &ranges[static_cast<std::size_t>(i)].second) != 2)
            return 1;
    }
    std::printf("%lld\n",
                static_cast<long long>(p0003::solve(n, q, ranges)));
    return 0;
}
