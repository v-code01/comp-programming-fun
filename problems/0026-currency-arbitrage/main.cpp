#include <cstddef>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 "n m"; then m lines "u v r" -- one unit of u buys r units of v.
// stdout: "YES" if an arbitrage loop exists, else "NO".
//
// rates arrive as reals with up to 6 decimals, so read them as doubles. m tops
// out at n*n = 90000, well inside scanf's reach.
int main() {
    int n = 0, m = 0;
    if (std::scanf("%d %d", &n, &m) != 2) return 1;

    std::vector<p0026::Offer> offers;
    offers.reserve(static_cast<std::size_t>(m));
    for (int e = 0; e < m; ++e) {
        int u = 0, v = 0;
        double r = 0.0;
        if (std::scanf("%d %d %lf", &u, &v, &r) != 3) return 1;
        offers.push_back({u, v, r});
    }

    std::printf("%s\n", p0026::has_arbitrage(n, offers) ? "YES" : "NO");
    return 0;
}
