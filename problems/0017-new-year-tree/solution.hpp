#pragma once

#include <bit>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0017 {

// New Year Tree. rooted tree, every vertex wears a color in 1..60. two ops:
//   repaint(v, c) -- paint every vertex in v's subtree the single color c.
//   distinct(v)   -- how many different colors live in v's subtree.
//
// THE FLATTEN. a subtree isn't a contiguous thing in a tree -- but an Euler tour
// makes it one. walk the tree depth-first, stamp each vertex with the time you
// first enter it (tin), and record the last entry-time seen anywhere below it
// (tout). then the subtree of v is exactly the vertices whose tin lands in
// [tin[v], tout[v]] -- one contiguous window on the tour. a subtree repaint is a
// range-assign; a subtree query is a range-read. the tree problem is now an array
// problem.
//
// WHY A 60-BIT MASK. colors cap at 60, and 60 bits ride inside one uint64. so
// represent the set of colors present over any window as a mask: bit c-1 set iff
// color c appears somewhere in there. combining two windows is a bitwise OR --
// union of sets, one instruction. "how many distinct colors" is then just
// popcount of the window's mask. no per-color bookkeeping, no sorting, no set
// container -- the whole distinct-count collapses to one 64-bit population count.
//
// THE SEGMENT TREE. over the tour positions 0..n-1. each node holds the OR-mask
// of its segment. range-assign to color c means "this whole window is now the
// single color c" -- overwrite the node's mask with the one-bit mask 1<<(c-1) and
// leave a lazy tag saying "everyone below me is this color too, push it when you
// descend." that lazy assignment is the load-bearing, bug-prone piece: a pending
// tag says the entire segment is uniformly one color, so pushing it down just
// copies the single-bit mask into both children and clears the parent's tag.
//
// COMPLEXITY. build is O(n). each repaint and each distinct touches O(log n)
// canonical nodes -- O((n + m) log n) over m ops.
//
// LOWER BOUND, honest. you must at least read the n colors, the n-1 edges, and
// the m queries -- Omega(n + m) just to consume the input, unavoidable. the
// segment tree adds a log n factor per op for the range-assign / range-union mix
// on arbitrary subtree windows; that log n is the standard price of dynamic range
// aggregation, not slack we could shave here.
class Tree {
 public:
    // color is 1-indexed: color[1..n] are the vertices' starting colors, color[0]
    // ignored. edges are 1-indexed undirected pairs. root is vertex 1.
    void build(int n, const std::vector<int>& color,
               const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        adj_.assign(static_cast<std::size_t>(n) + 1, {});
        for (const auto& e : edges) {
            adj_[static_cast<std::size_t>(e.first)].push_back(e.second);
            adj_[static_cast<std::size_t>(e.second)].push_back(e.first);
        }
        euler_tour(color);

        std::size_t cap = static_cast<std::size_t>(4 * (n_ > 0 ? n_ : 1));
        mask_.assign(cap, 0);
        lazy_.assign(cap, 0);
        if (n_ > 0) build_seg(1, 0, n_ - 1);
    }

    // paint v's subtree the single color c. range-assign over [tin, tout].
    void repaint(int v, int c) {
        std::uint64_t m = std::uint64_t{1} << (c - 1);  // one color, one bit.
        assign(1, 0, n_ - 1, tin_[static_cast<std::size_t>(v)],
               tout_[static_cast<std::size_t>(v)], m);
    }

    // number of distinct colors in v's subtree: popcount of the window's OR-mask.
    int distinct(int v) {
        std::uint64_t m = query(1, 0, n_ - 1, tin_[static_cast<std::size_t>(v)],
                                tout_[static_cast<std::size_t>(v)]);
        return std::popcount(m);
    }

 private:
    int n_ = 0;
    std::vector<std::vector<int>> adj_;
    std::vector<int> tin_;        // entry time per vertex -> subtree window start.
    std::vector<int> tout_;       // last entry time below -> window end.
    std::vector<int> euler_col_;  // color sitting at each tour position.
    std::vector<std::uint64_t> mask_;  // node OR-mask, segment-tree indexed.
    std::vector<std::uint64_t> lazy_;  // pending assignment: 0 = none, else the
                                       // single-bit mask the whole segment is set
                                       // to.

    // ITERATIVE Euler tour. a degenerate tree is a path of depth n; at n = 4e5 a
    // recursive DFS blows the call stack. so carry the recursion by hand: an
    // explicit vertex stack plus a per-vertex cursor into its adjacency list. we
    // enter a vertex, stamp tin, then advance its cursor one neighbor at a time;
    // when the cursor runs out we've finished the subtree, so stamp tout and pop.
    void euler_tour(const std::vector<int>& color) {
        std::size_t sz = static_cast<std::size_t>(n_) + 1;
        tin_.assign(sz, 0);
        tout_.assign(sz, 0);
        euler_col_.assign(static_cast<std::size_t>(n_), 0);
        std::vector<int> parent(sz, 0);
        std::vector<int> cursor(sz, 0);  // next neighbor index to inspect.
        std::vector<int> st;
        st.reserve(sz);

        int timer = 0;
        st.push_back(1);
        parent[1] = 0;
        tin_[1] = timer;
        euler_col_[static_cast<std::size_t>(timer)] = color[1];
        ++timer;

        while (!st.empty()) {
            int u = st.back();
            std::size_t us = static_cast<std::size_t>(u);
            if (cursor[us] < static_cast<int>(adj_[us].size())) {
                int w = adj_[us][static_cast<std::size_t>(cursor[us]++)];
                if (w == parent[us]) continue;  // don't walk back to the parent.
                parent[static_cast<std::size_t>(w)] = u;
                tin_[static_cast<std::size_t>(w)] = timer;
                euler_col_[static_cast<std::size_t>(timer)] =
                    color[static_cast<std::size_t>(w)];
                ++timer;
                st.push_back(w);
            } else {
                tout_[us] = timer - 1;  // everything below u is now stamped.
                st.pop_back();
            }
        }
    }

    void pull(int node) {
        mask_[static_cast<std::size_t>(node)] =
            mask_[static_cast<std::size_t>(2 * node)] |
            mask_[static_cast<std::size_t>(2 * node + 1)];
    }

    // stamp a whole segment as the single color m: its union is just m, and a
    // lazy tag records that every descendant is now m too.
    void apply(int node, std::uint64_t m) {
        mask_[static_cast<std::size_t>(node)] = m;
        lazy_[static_cast<std::size_t>(node)] = m;
    }

    void push_down(int node) {
        std::uint64_t m = lazy_[static_cast<std::size_t>(node)];
        if (m) {
            apply(2 * node, m);
            apply(2 * node + 1, m);
            lazy_[static_cast<std::size_t>(node)] = 0;
        }
    }

    void build_seg(int node, int lo, int hi) {
        if (lo == hi) {
            int c = euler_col_[static_cast<std::size_t>(lo)];
            mask_[static_cast<std::size_t>(node)] = std::uint64_t{1} << (c - 1);
            return;
        }
        int mid = (lo + hi) / 2;
        build_seg(2 * node, lo, mid);
        build_seg(2 * node + 1, mid + 1, hi);
        pull(node);
    }

    void assign(int node, int lo, int hi, int l, int r, std::uint64_t m) {
        if (r < lo || hi < l) return;            // segment outside the window.
        if (l <= lo && hi <= r) {                // fully inside: stamp and stop.
            apply(node, m);
            return;
        }
        push_down(node);
        int mid = (lo + hi) / 2;
        assign(2 * node, lo, mid, l, r, m);
        assign(2 * node + 1, mid + 1, hi, l, r, m);
        pull(node);
    }

    std::uint64_t query(int node, int lo, int hi, int l, int r) {
        if (r < lo || hi < l) return 0;  // disjoint contributes no colors.
        if (l <= lo && hi <= r) return mask_[static_cast<std::size_t>(node)];
        push_down(node);
        int mid = (lo + hi) / 2;
        return query(2 * node, lo, mid, l, r) |
               query(2 * node + 1, mid + 1, hi, l, r);
    }
};

}  // namespace p0017
