#pragma once

#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace p0057 {

// the dumb, obviously-correct oracle. no sidetracks, no path graph -- just grow
// every s->t walk-prefix in a min-heap and read off the lengths as they finish.
//
// a state is (length_so_far, node). seed (0, s). pop the smallest; if it sits on
// t, that prefix IS a completed s->t walk -- record its length. either way keep
// walking: push (length_so_far + w, x) for every out-edge (node,x,w), because a
// walk may pass through t and leave again. popping in non-decreasing length hands
// back walk lengths in exactly the order we want, multiplicity and ties honest.
//
// this is worst-case exponential and, on a graph with a zero-weight cycle, never
// terminates -- infinitely many length-tied walks. so it is for SMALL graphs
// only, and it carries a pop budget. `complete` says whether the answer is
// trustworthy: true if we collected k of them or the frontier emptied (fewer
// than k walks exist); false if we hit the budget first, in which case the
// caller must not compare against it.
class ReferenceSolver {
public:
    struct Result {
        std::vector<std::int64_t> lens;  // length k, padded with -1.
        bool complete;
    };

    // cap is the pop budget. tiny test graphs resolve k records long before it.
    Result solve(int n, int s, int t, int k, const std::vector<Edge>& edges,
                 long long cap = 5'000'000) {
        // forward adjacency -- plain vector of (head, weight) per vertex.
        std::vector<std::vector<std::pair<int, std::int64_t>>> adj(
            static_cast<std::size_t>(n) + 1);
        for (const auto& e : edges)
            adj[static_cast<std::size_t>(e.u)].push_back({e.v, e.w});

        using State = std::pair<std::int64_t, int>;  // (length_so_far, node)
        std::priority_queue<State, std::vector<State>, std::greater<>> pq;
        pq.push({0, s});

        std::vector<std::int64_t> rec;
        bool complete = false;
        long long pops = 0;
        while (!pq.empty()) {
            if (pops++ >= cap) break;  // budget blown -- answer is incomplete.
            const auto [len, node] = pq.top();
            pq.pop();
            if (node == t) {
                rec.push_back(len);
                if (static_cast<int>(rec.size()) == k) {
                    complete = true;
                    break;
                }
            }
            for (const auto& [x, w] : adj[static_cast<std::size_t>(node)])
                pq.push({len + w, x});
        }
        if (pq.empty()) complete = true;  // frontier drained: fewer than k walks.

        Result r;
        r.complete = complete;
        r.lens = rec;
        while (static_cast<int>(r.lens.size()) < k) r.lens.push_back(-1);
        return r;
    }
};

}  // namespace p0057
