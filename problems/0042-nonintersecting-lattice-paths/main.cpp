#include <cstdio>
#include <iostream>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 is k. then k lines "ax ay" (the starts), then k lines "bx by"
//         (the ends). the input is the compatible LGV setup -- both lists walk
//         south-east (x nondecreasing, y nonincreasing). see solution.hpp.
// stdout: the number of vertex-disjoint path families, mod 1e9+7.
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int k = 0;
    if (!(std::cin >> k)) return 1;
    if (k <= 0) {
        std::printf("0\n");
        return 0;
    }

    std::vector<p0042::Pt> a(static_cast<std::size_t>(k)), b(static_cast<std::size_t>(k));
    for (auto& p : a)
        if (!(std::cin >> p.x >> p.y)) return 1;
    for (auto& p : b)
        if (!(std::cin >> p.x >> p.y)) return 1;

    std::printf("%lld\n", static_cast<long long>(p0042::solve(a, b)));
    return 0;
}
