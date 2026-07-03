#pragma once

#include <algorithm>
#include <climits>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>

namespace p0020 {

// codeforces 863f, "almost permutation".
//
// systems framing: this is L2-norm load balancing. n items land in value-buckets
// 1..n, each item's bucket confined to a per-item interval [lo_i, hi_i] that the
// facts carve out. cost = sum_v cnt(v)^2 is the squared L2 norm of the load
// vector, and minimizing it minimizes load variance -- spread the n items across
// the buckets as evenly as the intervals allow. min sum of squares, subject to
// interval constraints, is min-cost flow.
//
// -----------------------------------------------------------------------------
// the convex-cost decomposition -- why 2k-1 edges reproduce cnt(v)^2
// -----------------------------------------------------------------------------
// sum of the first m odd numbers is m^2:  1 + 3 + 5 + ... + (2m-1) = m^2.
// (telescopes: k^2 - (k-1)^2 = 2k-1, sum k=1..m gives m^2 - 0.)
//
// so bucket v holding cnt(v) = m items costs 1 + 3 + ... + (2m-1). read that as a
// per-item ledger: the k-th item dropped into v pays exactly (2k-1). model it as
// n parallel unit-capacity edges val_v -> T with costs 1, 3, 5, ..., 2n-1. the
// cheapest unused edge on v is always the next odd number, so if MCMF ever routes
// m units through v it takes the m cheapest edges -- costs 1..(2m-1) -- summing to
// m^2. total flow cost = sum_v cnt(v)^2. exactly the objective, no gap.
//
// why the greedy successive-shortest-path is optimal here: the marginals 2k-1 are
// strictly increasing in k, so each bucket's cost is a convex function of its
// load. min-cost flow on a network whose only nontrivial costs are convex in
// per-arc flow is solved to optimality by SSP -- every augmentation pushes one
// unit down the globally cheapest residual path, and convexity means a unit never
// wants to displace an already-cheaper unit. no negative cycle survives. the odd
// edges *are* the convexity.
//
// -----------------------------------------------------------------------------
// the network
// -----------------------------------------------------------------------------
//   S -> pos_i            cap 1, cost 0        (each position ships one item)
//   pos_i -> val_v        cap 1, cost 0        for every v in [lo_i, hi_i]
//   val_v -> T            n unit edges, costs 1,3,5,...,2n-1   (convex bucket load)
// push n units S->T. answer = min cost. infeasible -> -1.
//
// feasibility: tighten lo/hi from the facts; if any lo_i > hi_i the interval is
// empty and no array exists -> -1. otherwise every position keeps >= 1 legal
// value and each bucket has capacity n >> n items, so a perfect S->T matching of
// all n units always exists (Hall is trivially satisfied). the flow < n guard
// below is therefore belt-and-suspenders, kept because the flow *is* the witness.
//
// -----------------------------------------------------------------------------
// complexity + the honest floor
// -----------------------------------------------------------------------------
// nodes: S + n positions + n values + T = 2n+2 = O(n).
// edges: n + O(n^2) + n^2 = O(n^2).
// units of flow: n. SSP runs n augmentations; each SPFA relaxation is O(V*E) =
// O(n * n^2) = O(n^3), so the solve is O(n^4). for n <= 50 that's ~6e6 relax
// steps -- sub-millisecond. it's polynomial, full stop; there is no linear floor
// it's pinned to. the only unavoidable Omega is reading the input: Omega(n + q)
// to see the array size and the q facts. the solve does real polynomial work on
// top of that, and honestly so -- MCMF, not a closed form.
//
// int width: costs sum to sum_v cnt(v)^2 <= (sum_v cnt(v))^2 = n^2 <= 2500 (all n
// items in one bucket is the max). path costs stay <= 2n-1 = 99. everything fits
// int with room to spare; the accumulator is int64 anyway, stated not assumed.
// -----------------------------------------------------------------------------

// one fact. 1-indexed, matching the statement. t==1: a_i >= v on [l,r].
// t==2: a_i <= v on [l,r].
struct Fact {
    int t, l, r, v;
};

// min-cost max-flow, successive shortest paths with SPFA. edges stored in pairs:
// forward at even index e, its residual partner at e^1. walking a path backward
// uses edges_[e^1].to for the tail node -- no separate parent array needed.
class MinCostFlow {
public:
    explicit MinCostFlow(int nodes) : head_(static_cast<std::size_t>(nodes), -1) {}

    void add_edge(int u, int v, int cap, int cost) {
        edges_.push_back({v, cap, cost, head_[static_cast<std::size_t>(u)]});
        head_[static_cast<std::size_t>(u)] = static_cast<int>(edges_.size()) - 1;
        edges_.push_back({u, 0, -cost, head_[static_cast<std::size_t>(v)]});
        head_[static_cast<std::size_t>(v)] = static_cast<int>(edges_.size()) - 1;
    }

    // push flow S->T until no augmenting path remains. returns {total flow, total
    // cost}. costs here are non-negative, so SPFA can't loop -- distances only
    // fall, and the residual graph has no negative cycle (convex marginals).
    std::pair<int, std::int64_t> run(int s, int t) {
        const int kInf = INT_MAX / 2;
        const int nodes = static_cast<int>(head_.size());
        int total_flow = 0;
        std::int64_t total_cost = 0;

        while (true) {
            std::vector<int> dist(static_cast<std::size_t>(nodes), kInf);
            std::vector<char> in_queue(static_cast<std::size_t>(nodes), 0);
            std::vector<int> pe(static_cast<std::size_t>(nodes), -1);  // path edge
            dist[static_cast<std::size_t>(s)] = 0;

            std::deque<int> q;
            q.push_back(s);
            in_queue[static_cast<std::size_t>(s)] = 1;
            while (!q.empty()) {
                int u = q.front();
                q.pop_front();
                in_queue[static_cast<std::size_t>(u)] = 0;
                for (int e = head_[static_cast<std::size_t>(u)]; e != -1;
                     e = edges_[static_cast<std::size_t>(e)].next) {
                    const Edge& ed = edges_[static_cast<std::size_t>(e)];
                    if (ed.cap <= 0) continue;
                    int nd = dist[static_cast<std::size_t>(u)] + ed.cost;
                    if (nd < dist[static_cast<std::size_t>(ed.to)]) {
                        dist[static_cast<std::size_t>(ed.to)] = nd;
                        pe[static_cast<std::size_t>(ed.to)] = e;
                        if (!in_queue[static_cast<std::size_t>(ed.to)]) {
                            in_queue[static_cast<std::size_t>(ed.to)] = 1;
                            q.push_back(ed.to);
                        }
                    }
                }
            }

            if (dist[static_cast<std::size_t>(t)] >= kInf) break;  // T unreachable.

            // bottleneck along the path, then apply it to both residuals.
            int push = INT_MAX;
            for (int v = t; v != s;) {
                int e = pe[static_cast<std::size_t>(v)];
                push = std::min(push, edges_[static_cast<std::size_t>(e)].cap);
                v = edges_[static_cast<std::size_t>(e ^ 1)].to;
            }
            for (int v = t; v != s;) {
                int e = pe[static_cast<std::size_t>(v)];
                edges_[static_cast<std::size_t>(e)].cap -= push;
                edges_[static_cast<std::size_t>(e ^ 1)].cap += push;
                v = edges_[static_cast<std::size_t>(e ^ 1)].to;
            }
            total_flow += push;
            total_cost += static_cast<std::int64_t>(push) *
                          dist[static_cast<std::size_t>(t)];
        }
        return {total_flow, total_cost};
    }

private:
    struct Edge {
        int to;
        int cap;
        int cost;
        int next;  // next edge out of the same tail, -1 to stop.
    };
    std::vector<Edge> edges_;
    std::vector<int> head_;
};

// the answer: minimum sum_v cnt(v)^2 over all arrays satisfying the facts, or -1.
inline std::int64_t solve(int n, const std::vector<Fact>& facts) {
    // tighten each position's interval. type-1 raises the floor, type-2 drops the
    // ceiling. start wide open at [1, n].
    std::vector<int> lo(static_cast<std::size_t>(n), 1);
    std::vector<int> hi(static_cast<std::size_t>(n), n);
    for (const Fact& f : facts) {
        for (int i = f.l - 1; i <= f.r - 1; ++i) {
            if (f.t == 1) {
                lo[static_cast<std::size_t>(i)] =
                    std::max(lo[static_cast<std::size_t>(i)], f.v);
            } else {
                hi[static_cast<std::size_t>(i)] =
                    std::min(hi[static_cast<std::size_t>(i)], f.v);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (lo[static_cast<std::size_t>(i)] > hi[static_cast<std::size_t>(i)]) {
            return -1;  // empty interval -- no array exists.
        }
    }

    // node layout: S=0, pos_i = 1+i, val_v = n+v (v in 1..n), T = 2n+1.
    const int kS = 0;
    const int kT = 2 * n + 1;
    MinCostFlow g(2 * n + 2);

    for (int i = 0; i < n; ++i) g.add_edge(kS, 1 + i, 1, 0);
    for (int i = 0; i < n; ++i) {
        for (int v = lo[static_cast<std::size_t>(i)];
             v <= hi[static_cast<std::size_t>(i)]; ++v) {
            g.add_edge(1 + i, n + v, 1, 0);
        }
    }
    // the convex bucket cost: k-th unit into value v pays (2k-1).
    for (int v = 1; v <= n; ++v) {
        for (int k = 1; k <= n; ++k) g.add_edge(n + v, kT, 1, 2 * k - 1);
    }

    auto [flow, cost] = g.run(kS, kT);
    if (flow < n) return -1;  // no full assignment -- infeasible (guard, see hdr).
    return cost;
}

}  // namespace p0020
