#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: line 1 is "n q". then q lines "t l r v". stdout: the min cost, or -1.
int main() {
    int n, q;
    if (std::scanf("%d %d", &n, &q) != 2) return 1;

    std::vector<p0020::Fact> facts(static_cast<std::size_t>(q));
    for (int i = 0; i < q; ++i) {
        p0020::Fact& f = facts[static_cast<std::size_t>(i)];
        if (std::scanf("%d %d %d %d", &f.t, &f.l, &f.r, &f.v) != 4) return 1;
    }

    std::printf("%lld\n",
                static_cast<long long>(p0020::solve(n, facts)));
    return 0;
}
