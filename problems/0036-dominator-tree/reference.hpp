#pragma once

#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

namespace p0036 {

// dumb and obviously correct, for tiny n only. it never mentions dfs numbers,
// semidominators, or link-eval -- it knows one fact and one fact only:
//
//   w dominates v  <=>  deleting w disconnects v from the root.
//
// so it deletes each vertex, re-runs a plain BFS, and watches who falls off the
// root's reachable set. that directly builds the dominator SET of every vertex.
// then idom(v) is the deepest proper dominator -- and because the dominators of
// v form a chain root ⊂ ... ⊂ idom(v), "deepest" is just "the one whose own
// dominator set is largest". no ties: a chain has a unique element per size.
//
// O(n * (n + m)) per query, exponential nowhere -- built to be believed, not to
// be fast. this is the ground truth the Lengauer-Tarjan output is diffed against.
class ReferenceSolver {
public:
    // returns idom[0..n], 1-indexed. idom[root] = 0, unreachable vertices 0.
    std::vector<int> solve(int n, const std::vector<std::pair<int, int>>& edges,
                           int root = 1) {
        n_ = n;
        // adjacency, self-loops dropped -- they never move reachability.
        adj_.assign(static_cast<std::size_t>(n) + 1, {});
        for (const auto& e : edges) {
            if (e.first == e.second) continue;
            if (e.first < 1 || e.first > n || e.second < 1 || e.second > n)
                continue;
            adj_[static_cast<std::size_t>(e.first)].push_back(e.second);
        }

        const std::vector<char> R = reach(root, -1);  // -1: remove nothing.

        // dom[v] as a bitmask over vertices 1..n (bit v-1). n stays tiny in
        // tests, so 64 bits is room to spare.
        std::vector<std::uint64_t> dom(static_cast<std::size_t>(n) + 1, 0);
        for (int v = 1; v <= n; ++v) {
            if (!R[static_cast<std::size_t>(v)]) continue;
            // root and v itself dominate v by definition; seed both.
            dom[static_cast<std::size_t>(v)] |= bit(root);
            dom[static_cast<std::size_t>(v)] |= bit(v);
        }
        // every other candidate d: remove it, see who the root loses.
        for (int d = 1; d <= n; ++d) {
            if (d == root) continue;
            if (!R[static_cast<std::size_t>(d)]) continue;  // unreachable d dominates nobody.
            const std::vector<char> Rd = reach(root, d);
            for (int v = 1; v <= n; ++v) {
                if (v == d) continue;
                if (!R[static_cast<std::size_t>(v)]) continue;
                if (!Rd[static_cast<std::size_t>(v)])
                    dom[static_cast<std::size_t>(v)] |= bit(d);  // d dominates v.
            }
        }

        std::vector<int> idom(static_cast<std::size_t>(n) + 1, 0);
        for (int v = 1; v <= n; ++v) {
            if (v == root || !R[static_cast<std::size_t>(v)]) continue;
            // among proper dominators, pick the one with the biggest dom set --
            // the deepest link in the chain.
            int best = 0, best_size = -1;
            std::uint64_t mask = dom[static_cast<std::size_t>(v)] & ~bit(v);
            while (mask) {
                const int d = lsb_index(mask) + 1;
                mask &= mask - 1;
                const int sz = popcount(dom[static_cast<std::size_t>(d)]);
                if (sz > best_size) {
                    best_size = sz;
                    best = d;
                }
            }
            idom[static_cast<std::size_t>(v)] = best;
        }
        return idom;
    }

private:
    static std::uint64_t bit(int v) {
        return std::uint64_t{1} << (v - 1);
    }
    static int popcount(std::uint64_t x) {
        return __builtin_popcountll(x);
    }
    static int lsb_index(std::uint64_t x) {
        return __builtin_ctzll(x);
    }

    // BFS from root over the graph with `removed` deleted (removed = -1 keeps
    // everything). returns a 1-indexed reached flag.
    std::vector<char> reach(int root, int removed) {
        std::vector<char> seen(static_cast<std::size_t>(n_) + 1, 0);
        if (root == removed) return seen;  // root gone -> nothing reachable.
        std::queue<int> q;
        seen[static_cast<std::size_t>(root)] = 1;
        q.push(root);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            for (const int w : adj_[static_cast<std::size_t>(u)]) {
                if (w == removed) continue;
                if (!seen[static_cast<std::size_t>(w)]) {
                    seen[static_cast<std::size_t>(w)] = 1;
                    q.push(w);
                }
            }
        }
        return seen;
    }

    int n_ = 0;
    std::vector<std::vector<int>> adj_;
};

}  // namespace p0036
