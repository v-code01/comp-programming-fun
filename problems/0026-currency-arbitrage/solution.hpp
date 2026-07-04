#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace p0026 {

// one exchange offer -- a single unit of u buys `rate` units of v. rate > 0.
struct Offer {
    int u;
    int v;
    double rate;
};

// currency arbitrage. is there a loop of conversions that ends holding strictly
// more than it started? triangular / statistical arb -- turn a dollar into euros
// into pounds into dollars and walk away with $1.03, risk-free, if the rates line
// up wrong.
//
// a loop pays when the product of its rates clears 1. take -log of every rate and
// that product collapses into a sum:
//   prod(rate) > 1  <=>  sum(log rate) > 0  <=>  sum(-log rate) < 0
// so weight w(u->v) = -log(rate) rewrites "arbitrage cycle" as "negative-weight
// cycle" -- the same object, exactly. find one and you've found the arb.
//
// bellman-ford does the finding. seed every node at dist 0. that's a virtual
// super-source wired to each node with a free edge, so every currency is a
// candidate start at once. relax all m edges n-1 times -- enough to settle every
// shortest path a negative cycle can't reach. then sweep once more: an edge that
// still relaxes has no lower bound on its distance, and that only happens with a
// negative cycle upstream. that edge is the arbitrage.
//
// complexity: O(n * m). the lower bound is honest either way -- reading the graph
// is Omega(n + m), and the n-1 sweeps over m edges are the work stacked on top.
// eps guards the seam. a loop whose product sits a float-hair off 1 must not flip
// the verdict on rounding noise, so a change under eps doesn't count as a change.
inline bool has_arbitrage(int n, const std::vector<Offer>& offers) {
    constexpr double kEps = 1e-9;

    // dist seeded at 0 for every node -- the super-source trick, no extra vertex.
    std::vector<double> dist(static_cast<std::size_t>(n) + 1, 0.0);

    // -log once, up front. log is the expensive call here, not the compares.
    std::vector<double> w;
    w.reserve(offers.size());
    for (const Offer& o : offers) w.push_back(-std::log(o.rate));

    // n-1 relaxation sweeps. bail the instant a sweep moves nothing by more than
    // eps -- the distances have settled and the detector can only agree. a real
    // negative cycle keeps paying out every sweep, so it never triggers this.
    for (int pass = 0; pass < n - 1; ++pass) {
        bool changed = false;
        for (std::size_t e = 0; e < offers.size(); ++e) {
            const std::size_t u = static_cast<std::size_t>(offers[e].u);
            const std::size_t v = static_cast<std::size_t>(offers[e].v);
            const double cand = dist[u] + w[e];
            if (cand < dist[v] - kEps) {
                dist[v] = cand;
                changed = true;
            }
        }
        if (!changed) return false;
    }

    // the n-th sweep is the detector. one more relaxable edge -- negative cycle.
    for (std::size_t e = 0; e < offers.size(); ++e) {
        const std::size_t u = static_cast<std::size_t>(offers[e].u);
        const std::size_t v = static_cast<std::size_t>(offers[e].v);
        if (dist[u] + w[e] < dist[v] - kEps) return true;
    }
    return false;
}

}  // namespace p0026
