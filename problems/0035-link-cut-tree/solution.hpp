#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0035 {

// the link-cut tree -- the hardest standard data structure, and the reason it's
// here. a forest, fully dynamic: link two trees, cut an edge, add a constant to
// every vertex on a path, sum a path. all four in O(log n) amortized.
//
// THE IDEA. decompose the represented forest into preferred paths. each preferred
// path is stored as a splay tree keyed by depth -- shallow on the left, deep on
// the right. a splay tree's in-order traversal is exactly the path from its
// top vertex down to its bottom vertex. the splay trees are stitched together by
// path-parent pointers: the root of one splay tree points up to a vertex in the
// splay tree above it, but is NOT that vertex's child. that asymmetry is the
// whole trick, and it's what isRoot() below tests.
//
// access(x) -- expose the root-to-x path. walk from x up to the true root,
// splaying at each step and re-hanging the preferred child, so afterward one
// splay tree holds exactly root..x and x sits at its top. O(log n) amortized by
// the same heavy-path / access-lemma potential that pays for splay itself.
//
// makeRoot(x) -- reroot the represented tree at x. access(x) puts root..x in one
// splay tree; a lazy REVERSE flag on that tree flips depth order, so x becomes the
// shallowest vertex -- the new root. evert is nothing but a reverse.
//
// with makeRoot, every path op reduces to "put u..v alone in one splay tree":
//   pathSum(u,v)   = makeRoot(u); access(v); read sum at the splay root.
//   pathAdd(u,v,x) = makeRoot(u); access(v); lazy-add x to the whole splay tree.
//   link(u,v)      = makeRoot(u); point u's path-parent at v.
//   cut(u,v)       = makeRoot(u); access(v); u is v's shallower child -- detach.
//
// THE BUG-PRONE SPOT -- push_down carrying BOTH reverse and add. each node holds
// a pending reverse flag and a pending add. push_down ships both to the children.
// they commute: reverse only swaps a node's two children (re-orders the path),
// while add is uniform across the whole subtree -- adding then reversing lands the
// same values as reversing then adding, and the subtree SUM is invariant under
// reversal. so push_down applies them in either order and the aggregates stay
// exact through every rotation. get that composition wrong and makeRoot silently
// corrupts a path add; the brute oracle in reference.hpp is the gate that catches
// it.
//
// COMPLEXITY. O((n + q) log n) amortized total. reading input alone is
// Omega(n + q), so the log n per op is the only slack -- and the naive
// alternative, walking the path vertex by vertex, is O(n) per op. on a chain of
// 1e5 with 1e5 path queries that's 1e10 vs the LCT's ~1e5 * 17. the LCT wins by
// the depth of the tree.
//
// INT64, and why the sums fit. a vertex value is a_i plus every x from a path-add
// that covered it: |a_i|, |x| <= 1e9. a path holds at most L+1 vertices where L is
// the number of links so far, and each of A add-ops touches at most L+1 vertices.
// links and adds are both ops, so A + L <= q = 1e5. a query sum is bounded by
//   |base| + |adds| <= (L+1)*1e9 + A*(L+1)*1e9 <= 1e14 + 1e9 * A*(L+1).
// A*(L+1) peaks near (5e4)^2 = 2.5e9 under A + L <= 1e5, so the sum stays under
// ~2.5e18 -- inside int64's 9.22e18. the same budget bounds add*size during a
// push_down, since add[t]*sz[t] is just the slice of the true sum being carried.
class LinkCutTree {
public:
    // build over a_i for vertices 1..n. node 0 is the null sentinel -- size 0,
    // sum 0 -- so child reads off the bottom of a splay tree fold in as nothing.
    void build(const std::vector<std::int64_t>& a) {
        n_ = static_cast<int>(a.size());
        std::size_t cap = static_cast<std::size_t>(n_ + 1);
        ch_.assign(cap, {0, 0});
        fa_.assign(cap, 0);
        val_.assign(cap, 0);
        sum_.assign(cap, 0);
        add_.assign(cap, 0);
        sz_.assign(cap, 0);
        rev_.assign(cap, 0);
        for (int i = 1; i <= n_; ++i) {
            val_[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i - 1)];
            sum_[static_cast<std::size_t>(i)] = val_[static_cast<std::size_t>(i)];
            sz_[static_cast<std::size_t>(i)] = 1;
        }
        stk_.reserve(cap);
    }

    // add an edge u--v. u and v live in different trees by contract, so rerooting
    // u and pointing its path-parent at v joins the two without a cycle.
    void link(int u, int v) {
        makeRoot(u);
        fa_[static_cast<std::size_t>(u)] = v;  // path-parent, not a splay child.
    }

    // remove the edge u--v. it exists by contract, so after makeRoot(u) and
    // access(v) the whole path u..v is the two-node splay tree {u, v} with u as v's
    // shallower (left) child -- snip that link.
    void cut(int u, int v) {
        makeRoot(u);
        access(v);
        // v is the splay root; u is its left child and its only other node.
        ch_[static_cast<std::size_t>(v)][0] = 0;
        fa_[static_cast<std::size_t>(u)] = 0;
        pull(v);
    }

    // x += ... on every vertex of the path u..v. makeRoot + access isolate the
    // path as v's splay subtree, then one lazy add covers all of it.
    void pathAdd(int u, int v, std::int64_t x) {
        makeRoot(u);
        access(v);
        apply_add(v, x);
    }

    // sum of values over the path u..v. same isolation, then read the fold.
    std::int64_t pathSum(int u, int v) {
        makeRoot(u);
        access(v);
        return sum_[static_cast<std::size_t>(v)];
    }

private:
    int n_ = 0;
    std::vector<std::array<int, 2>> ch_;  // splay children: [0]=shallower, [1]=deeper.
    std::vector<int> fa_;                 // splay parent OR path-parent (see isRoot).
    std::vector<std::int64_t> val_;       // the vertex's own value.
    std::vector<std::int64_t> sum_;       // sum over this splay subtree.
    std::vector<std::int64_t> add_;       // pending add for the subtree.
    std::vector<int> sz_;                 // subtree size -- turns add into add*sz on sum.
    std::vector<char> rev_;               // pending reverse for the subtree.
    std::vector<int> stk_;                // scratch for the top-down push before a splay.

    // x is the top of its splay tree iff its parent doesn't own it as a child.
    // a path-parent points up but the target's children don't point back -- that's
    // the boundary between splay trees, and the definition of a preferred path's root.
    bool isRoot(int x) const {
        int f = fa_[static_cast<std::size_t>(x)];
        return ch_[static_cast<std::size_t>(f)][0] != x &&
               ch_[static_cast<std::size_t>(f)][1] != x;
    }

    // recombine a node from its children. size and sum only -- val is the node's own.
    void pull(int x) {
        int l = ch_[static_cast<std::size_t>(x)][0];
        int r = ch_[static_cast<std::size_t>(x)][1];
        sz_[static_cast<std::size_t>(x)] =
            sz_[static_cast<std::size_t>(l)] + sz_[static_cast<std::size_t>(r)] + 1;
        sum_[static_cast<std::size_t>(x)] = sum_[static_cast<std::size_t>(l)] +
                                            sum_[static_cast<std::size_t>(r)] +
                                            val_[static_cast<std::size_t>(x)];
    }

    // stage a reverse on x's subtree: swap its children now, remember to push the
    // flag down later. sum is reversal-invariant, so it doesn't move.
    void apply_rev(int x) {
        if (x == 0) return;  // never toggle the sentinel.
        std::swap(ch_[static_cast<std::size_t>(x)][0],
                  ch_[static_cast<std::size_t>(x)][1]);
        rev_[static_cast<std::size_t>(x)] ^= 1;
    }

    // stage an add on x's subtree: every value shifts by v, so the fold shifts by
    // v*size, and the pending add accumulates.
    void apply_add(int x, std::int64_t v) {
        if (x == 0) return;  // guarding this keeps the sentinel's val at zero.
        val_[static_cast<std::size_t>(x)] += v;
        sum_[static_cast<std::size_t>(x)] +=
            v * static_cast<std::int64_t>(sz_[static_cast<std::size_t>(x)]);
        add_[static_cast<std::size_t>(x)] += v;
    }

    // ship both pending flags to the children. reverse and add commute here, so
    // the order between these two blocks is free -- this is the composition the
    // whole structure hinges on.
    void push(int x) {
        std::size_t xs = static_cast<std::size_t>(x);
        if (rev_[xs]) {
            apply_rev(ch_[xs][0]);
            apply_rev(ch_[xs][1]);
            rev_[xs] = 0;
        }
        if (add_[xs] != 0) {
            apply_add(ch_[xs][0], add_[xs]);
            apply_add(ch_[xs][1], add_[xs]);
            add_[xs] = 0;
        }
    }

    // one splay rotation, path-parent aware. if y was a splay root, z is its
    // path-parent: x inherits that pointer, but z must NOT gain x as a child --
    // the isRoot(y) guard is exactly that case.
    void rotate(int x) {
        int y = fa_[static_cast<std::size_t>(x)];
        int z = fa_[static_cast<std::size_t>(y)];
        int k = (ch_[static_cast<std::size_t>(y)][1] == x) ? 1 : 0;  // x is y's k-child.
        if (!isRoot(y)) {
            int side = (ch_[static_cast<std::size_t>(z)][1] == y) ? 1 : 0;
            ch_[static_cast<std::size_t>(z)][side] = x;
        }
        fa_[static_cast<std::size_t>(x)] = z;
        int b = ch_[static_cast<std::size_t>(x)][k ^ 1];  // x's inner child moves to y.
        ch_[static_cast<std::size_t>(y)][k] = b;
        if (b) fa_[static_cast<std::size_t>(b)] = y;
        ch_[static_cast<std::size_t>(x)][k ^ 1] = y;
        fa_[static_cast<std::size_t>(y)] = x;
        pull(y);  // y sank -- refold it first.
        pull(x);  // then x above it.
    }

    // splay x to the top of its splay tree. first push every pending flag from the
    // splay root down to x -- rotations read children, so the pending flags have
    // to be resolved along the path before a single rotate fires. then the usual
    // zig-zig / zig-zag climb.
    void splay(int x) {
        stk_.clear();
        int y = x;
        stk_.push_back(y);
        while (!isRoot(y)) {
            y = fa_[static_cast<std::size_t>(y)];
            stk_.push_back(y);
        }
        for (int i = static_cast<int>(stk_.size()) - 1; i >= 0; --i)
            push(stk_[static_cast<std::size_t>(i)]);

        while (!isRoot(x)) {
            y = fa_[static_cast<std::size_t>(x)];
            int z = fa_[static_cast<std::size_t>(y)];
            if (!isRoot(y)) {
                // same side -> rotate the middle first (zig-zig); else rotate x (zig-zag).
                bool same = (ch_[static_cast<std::size_t>(z)][1] == y) ==
                            (ch_[static_cast<std::size_t>(y)][1] == x);
                rotate(same ? y : x);
            }
            rotate(x);
        }
    }

    // expose the path from the represented root down to x, then splay x to the top
    // of the resulting splay tree. each iteration splays the current vertex, swaps
    // in the path we came from as its deeper (preferred) child, refolds, and climbs
    // the path-parent pointer.
    void access(int x) {
        int last = 0;
        for (int y = x; y != 0; y = fa_[static_cast<std::size_t>(y)]) {
            splay(y);
            ch_[static_cast<std::size_t>(y)][1] = last;  // detach old preferred, attach lower.
            pull(y);
            last = y;
        }
        splay(x);
    }

    // reroot at x: expose root..x, then reverse it so x -- formerly the deepest --
    // becomes the shallowest, i.e. the new root.
    void makeRoot(int x) {
        access(x);
        apply_rev(x);
    }
};

}  // namespace p0035
