#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0019 {

// the dumb, obviously-correct oracle -- exponential, tiny n only. it shares not a
// line of logic with the greedy: no tails, no peaks, no sqrt trick. it enumerates
// every length-k path and then brute-forces the maximum set of pairwise
// vertex-disjoint ones. that independence is the whole point -- it is the ground
// truth the greedy has to match on every tree and every k.
//
// two pieces:
//   1. every simple path of exactly k vertices. in a tree a simple path is fixed
//      by its two endpoints, so a path is just a pair (u, w) whose unique route
//      spans k vertices. we record each as a bitmask of its vertices.
//   2. maximum vertex-disjoint packing of those masks. classic exact exponential
//      set packing: look at the lowest free vertex and either leave it out of
//      every path or cover it with one path through it, recurse, memoize on the
//      used-set. every vertex is covered-or-not, so this misses no packing.
//
// n must stay small -- the memo is indexed by a 2^n bitmask. the differential test
// feeds it trees up to ~12 vertices, where it answers instantly.
class ReferenceSolver {
public:
    // f(k) for every k in [1, n], as ans[k-1]. matches the fast solver's shape.
    std::vector<int> solve(int n,
                           const std::vector<std::pair<int, int>>& edges) {
        std::vector<int> ans(static_cast<std::size_t>(n < 0 ? 0 : n), 0);
        if (n <= 0) return ans;
        n_ = n;
        const std::size_t N = static_cast<std::size_t>(n);

        // plain adjacency lists -- clarity over speed, the trees are tiny.
        std::vector<std::vector<int>> g(N);
        for (const auto& e : edges) {
            g[static_cast<std::size_t>(e.first - 1)].push_back(e.second - 1);
            g[static_cast<std::size_t>(e.second - 1)].push_back(e.first - 1);
        }

        // all-pairs distance and parent, one BFS per source. O(n^2) on a tiny n.
        std::vector<std::vector<int>> dist(N, std::vector<int>(N, -1));
        std::vector<std::vector<int>> par(N, std::vector<int>(N, -1));
        for (std::size_t s = 0; s < N; ++s) {
            std::vector<int> q;
            q.reserve(N);
            q.push_back(static_cast<int>(s));
            dist[s][s] = 0;
            for (std::size_t qi = 0; qi < q.size(); ++qi) {
                const int u = q[qi];
                for (const int w : g[static_cast<std::size_t>(u)]) {
                    if (dist[s][static_cast<std::size_t>(w)] == -1) {
                        dist[s][static_cast<std::size_t>(w)] =
                            dist[s][static_cast<std::size_t>(u)] + 1;
                        par[s][static_cast<std::size_t>(w)] = u;
                        q.push_back(w);
                    }
                }
            }
        }

        for (int k = 1; k <= n; ++k) {
            // gather the length-k path masks.
            std::vector<std::uint32_t> paths;
            if (k == 1) {
                for (int v = 0; v < n; ++v)
                    paths.push_back(std::uint32_t{1} << v);
            } else {
                for (int u = 0; u < n; ++u)
                    for (int w = u + 1; w < n; ++w) {
                        if (dist[static_cast<std::size_t>(u)]
                                [static_cast<std::size_t>(w)] != k - 1)
                            continue;
                        // walk w -> u along parents, collecting the vertices.
                        std::uint32_t mask = 0;
                        int x = w;
                        while (x != -1) {
                            mask |= std::uint32_t{1} << x;
                            if (x == u) break;
                            x = par[static_cast<std::size_t>(u)]
                                   [static_cast<std::size_t>(x)];
                        }
                        paths.push_back(mask);
                    }
            }

            // maximum disjoint packing over those masks.
            memo_.assign(std::size_t{1} << N, -1);
            ans[static_cast<std::size_t>(k - 1)] = pack(0, paths);
        }
        return ans;
    }

private:
    // max disjoint paths using only vertices outside `used`. decide the lowest
    // free vertex: skip it, or cover it with one path through it. memoized on
    // `used` -- the answer depends on nothing else.
    int pack(std::uint32_t used, const std::vector<std::uint32_t>& paths) {
        const int cached = memo_[used];
        if (cached != -1) return cached;

        int v = -1;
        for (int i = 0; i < n_; ++i)
            if (!((used >> i) & 1u)) {
                v = i;
                break;
            }
        if (v < 0) return memo_[used] = 0;  // everything placed.

        // leave v out of every path.
        int best = pack(used | (std::uint32_t{1} << v), paths);
        // or cover v with a path through it that avoids the used set.
        for (const std::uint32_t p : paths) {
            if (((p >> v) & 1u) && !(p & used)) {
                const int cand = 1 + pack(used | p, paths);
                if (cand > best) best = cand;
            }
        }
        return memo_[used] = best;
    }

    int n_ = 0;
    std::vector<int> memo_;  // memo_[used] = best packing, sized 2^n.
};

// free-function form, for one-shot calls.
inline std::vector<int> you_are_given_a_tree_ref(
    int n, const std::vector<std::pair<int, int>>& edges) {
    ReferenceSolver ref;
    return ref.solve(n, edges);
}

}  // namespace p0019
