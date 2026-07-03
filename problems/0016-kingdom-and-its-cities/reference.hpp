#pragma once

#include <utility>
#include <vector>

namespace p0016 {

// dumb and obviously correct, for tiny trees only. it does not compress anything
// or reason about junctions -- it just tries every possible set of non-important
// cities to remove and keeps the smallest one that works.
//
//   -1: if any tree edge joins two important cities, no removal can split them
//       (you may only remove non-important cities). return -1.
//   else: over all 2^(n-k) subsets of the non-important cities, remove that
//       subset and check -- by union-find on the surviving edges -- that no two
//       important cities share a component. the answer is the smallest subset
//       size that passes.
//
// exponential in the number of non-important cities, so keep n tiny. that's the
// whole point: it's the ground truth the virtual-tree dp is diffed against.
class ReferenceSolver {
public:
    // nodes 1..n, edges 0-indexed pairs, important = the query's important cities.
    int solve(int n, const std::vector<std::pair<int, int>>& edges,
              const std::vector<int>& important) {
        std::vector<char> is_imp(static_cast<std::size_t>(n) + 1, 0);
        for (int a : important) is_imp[static_cast<std::size_t>(a)] = 1;

        // two adjacent important cities -> impossible.
        for (const auto& e : edges) {
            if (is_imp[static_cast<std::size_t>(e.first)] &&
                is_imp[static_cast<std::size_t>(e.second)]) {
                return -1;
            }
        }

        // the cities we're allowed to remove.
        std::vector<int> free;
        for (int v = 1; v <= n; ++v) {
            if (!is_imp[static_cast<std::size_t>(v)]) free.push_back(v);
        }
        const int m = static_cast<int>(free.size());

        int best = n + 1;  // removing everything free always works, so this bounds it.
        for (int mask = 0; mask < (1 << m); ++mask) {
            const int bits = popcount(mask);
            if (bits >= best) continue;  // can't beat the current best -- skip.

            // mark this removed set.
            std::vector<char> removed(static_cast<std::size_t>(n) + 1, 0);
            for (int i = 0; i < m; ++i) {
                if (mask & (1 << i)) removed[static_cast<std::size_t>(free[static_cast<std::size_t>(i)])] = 1;
            }

            // union the surviving graph, then check important cities are all in
            // distinct components.
            parent_.assign(static_cast<std::size_t>(n) + 1, 0);
            for (int v = 0; v <= n; ++v) parent_[static_cast<std::size_t>(v)] = v;
            for (const auto& e : edges) {
                if (removed[static_cast<std::size_t>(e.first)] ||
                    removed[static_cast<std::size_t>(e.second)]) {
                    continue;
                }
                unite(e.first, e.second);
            }

            bool ok = true;
            std::vector<char> comp_seen(static_cast<std::size_t>(n) + 1, 0);
            for (int a : important) {
                const int r = find(a);
                if (comp_seen[static_cast<std::size_t>(r)]) {
                    ok = false;  // two important cities in one component -- still joined.
                    break;
                }
                comp_seen[static_cast<std::size_t>(r)] = 1;
            }
            if (ok) best = bits;
        }
        return best;
    }

private:
    static int popcount(int x) {
        int c = 0;
        while (x) {
            x &= x - 1;
            ++c;
        }
        return c;
    }

    int find(int x) {
        while (parent_[static_cast<std::size_t>(x)] != x) {
            parent_[static_cast<std::size_t>(x)] =
                parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(x)])];
            x = parent_[static_cast<std::size_t>(x)];
        }
        return x;
    }

    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent_[static_cast<std::size_t>(a)] = b;
    }

    std::vector<int> parent_;
};

// free-function convenience matching the solution's shape.
inline int reference(int n, const std::vector<std::pair<int, int>>& edges,
                     const std::vector<int>& important) {
    ReferenceSolver s;
    return s.solve(n, edges, important);
}

}  // namespace p0016
