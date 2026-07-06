#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0035 {

// the oracle. an adjacency list and a value array, nothing else. link and cut
// edit the adjacency; a path op walks the path vertex by vertex and loops it.
// O(n) per op, obviously correct, no splay anything -- that's the point. the LCT
// has to agree with this on every input. same API (build/link/cut/pathAdd/pathSum)
// so one op program replays on both and the type-4 streams get diffed.
class Reference {
public:
    // vertices are 1..n; slot 0 is unused so indices match the LCT's 1-based nodes.
    void build(const std::vector<std::int64_t>& a) {
        n_ = static_cast<int>(a.size());
        val_.assign(static_cast<std::size_t>(n_ + 1), 0);
        for (int i = 1; i <= n_; ++i)
            val_[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i - 1)];
        adj_.assign(static_cast<std::size_t>(n_ + 1), {});
    }

    void link(int u, int v) {
        adj_[static_cast<std::size_t>(u)].push_back(v);
        adj_[static_cast<std::size_t>(v)].push_back(u);
    }

    // drop the edge from both endpoints' lists. the edge exists by contract.
    void cut(int u, int v) {
        erase(u, v);
        erase(v, u);
    }

    void pathAdd(int u, int v, std::int64_t x) {
        for (int c : path(u, v)) val_[static_cast<std::size_t>(c)] += x;
    }

    std::int64_t pathSum(int u, int v) {
        std::int64_t s = 0;
        for (int c : path(u, v)) s += val_[static_cast<std::size_t>(c)];
        return s;
    }

private:
    int n_ = 0;
    std::vector<std::int64_t> val_;
    std::vector<std::vector<int>> adj_;

    void erase(int a, int b) {
        auto& L = adj_[static_cast<std::size_t>(a)];
        for (std::size_t i = 0; i < L.size(); ++i)
            if (L[i] == b) {
                L.erase(L.begin() + static_cast<std::ptrdiff_t>(i));
                return;
            }
    }

    // BFS from u, stop at v, then read parents back to recover the u..v path.
    // trees are acyclic and u,v are connected by contract, so the path is unique.
    std::vector<int> path(int u, int v) {
        std::vector<int> par(static_cast<std::size_t>(n_ + 1), 0);
        std::vector<char> seen(static_cast<std::size_t>(n_ + 1), 0);
        std::vector<int> q;
        q.push_back(u);
        seen[static_cast<std::size_t>(u)] = 1;
        for (std::size_t i = 0; i < q.size(); ++i) {
            int x = q[i];
            if (x == v) break;
            for (int w : adj_[static_cast<std::size_t>(x)]) {
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                par[static_cast<std::size_t>(w)] = x;
                q.push_back(w);
            }
        }
        // par[u] stayed 0 -- the sentinel that ends the climb. handles u == v too.
        std::vector<int> p;
        for (int c = v; c != 0; c = par[static_cast<std::size_t>(c)]) p.push_back(c);
        return p;
    }
};

}  // namespace p0035
