#pragma once

#include <bit>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0017 {

// the oracle. same Euler-tour flatten, but no cleverness on top: a plain array of
// per-position colors and the definitions read straight off the statement.
//   repaint(v, c) -- loop v's window and write c into every cell. O(subtree).
//   distinct(v)   -- loop v's window, OR each color's bit into a mask, popcount.
// O(n) per op, obviously correct, nothing to get wrong. the segment tree has to
// agree with this on every input -- that's the whole reason it exists.
class Reference {
 public:
    void build(int n, const std::vector<int>& color,
               const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        std::vector<std::vector<int>> adj(static_cast<std::size_t>(n) + 1);
        for (const auto& e : edges) {
            adj[static_cast<std::size_t>(e.first)].push_back(e.second);
            adj[static_cast<std::size_t>(e.second)].push_back(e.first);
        }

        std::size_t sz = static_cast<std::size_t>(n) + 1;
        tin_.assign(sz, 0);
        tout_.assign(sz, 0);
        col_.assign(static_cast<std::size_t>(n), 0);
        std::vector<int> parent(sz, 0);
        std::vector<int> cursor(sz, 0);
        std::vector<int> st;
        st.reserve(sz);

        int timer = 0;
        st.push_back(1);
        parent[1] = 0;
        tin_[1] = timer;
        col_[static_cast<std::size_t>(timer)] = color[1];
        ++timer;

        while (!st.empty()) {
            int u = st.back();
            std::size_t us = static_cast<std::size_t>(u);
            if (cursor[us] < static_cast<int>(adj[us].size())) {
                int w = adj[us][static_cast<std::size_t>(cursor[us]++)];
                if (w == parent[us]) continue;
                parent[static_cast<std::size_t>(w)] = u;
                tin_[static_cast<std::size_t>(w)] = timer;
                col_[static_cast<std::size_t>(timer)] =
                    color[static_cast<std::size_t>(w)];
                ++timer;
                st.push_back(w);
            } else {
                tout_[us] = timer - 1;
                st.pop_back();
            }
        }
    }

    void repaint(int v, int c) {
        for (int p = tin_[static_cast<std::size_t>(v)];
             p <= tout_[static_cast<std::size_t>(v)]; ++p) {
            col_[static_cast<std::size_t>(p)] = c;
        }
    }

    int distinct(int v) {
        std::uint64_t m = 0;
        for (int p = tin_[static_cast<std::size_t>(v)];
             p <= tout_[static_cast<std::size_t>(v)]; ++p) {
            m |= std::uint64_t{1} << (col_[static_cast<std::size_t>(p)] - 1);
        }
        return std::popcount(m);
    }

 private:
    int n_ = 0;
    std::vector<int> tin_;
    std::vector<int> tout_;
    std::vector<int> col_;  // current color at each tour position.
};

}  // namespace p0017
