#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0013 {

// rooted tree at vertex 1, n vertices, vertex v painted color c_v. in a subtree
// a color "dominates" when nothing outnumbers it -- ties allowed. for every v,
// sum the colors that hit the max count inside subtree(v). (codeforces 600E,
// n <= 1e5, c_v in [1, n], so a subtree answer can reach ~n distinct colors of
// size ~n -- 1e10, an int64.)
//
// the naive answer recomputes each subtree from scratch: O(n) colors per vertex,
// n vertices, O(n^2). the whole subtree of one vertex overlaps enormously with
// its parent's -- recomputing throws that away. dsu on tree (small-to-large,
// "sack") keeps the overlap and pays for the rest.
//
// THE HEAVY-CHILD SPLIT. root the tree. call the child with the largest subtree
// the heavy child; the rest are light. the heavy edge of each vertex marks one
// path down; every other edge is light. walk from any vertex up to the root and
// count the light edges you cross -- each light edge steps you into a subtree at
// most half the size of the one above (its sibling heavy subtree is at least as
// big, so the parent is >= 2x the light child). size halves at most n times
// before it hits one, so any root-to-vertex path crosses O(log n) light edges.
//
// THE SACK RECURSION. process a vertex's light children first, each time wiping
// the color counts clean when it returns. then process the heavy child and KEEP
// its counts live. now cnt[] already holds the heavy subtree for free -- fold in
// the vertex itself and each light subtree by brute, read off the answer, and if
// this vertex is itself light, wipe everything so the sibling starts clean.
//
// THE CHARGE. a vertex's colors get re-inserted once for every light edge on its
// path to the root -- that's the "fold the light subtree back in" step, which
// touches every descendant of a light child. O(log n) light edges above it means
// O(log n) insertions total, summed over all vertices O(n log n). the heavy
// subtree is never re-inserted; that's where the saving lives.
//
// LOWER BOUND. reading the tree is Omega(n) -- n colors and n-1 edges, each must
// be looked at once. the log factor on top is the small-to-large charge above:
// a balanced binary tree forces every vertex through log n light edges, so
// Theta(n log n) insertions happen and the bound is tight, not loose.
//
// EULER FLATTENING makes "fold in a whole subtree" a contiguous scan instead of
// a recursion. a dfs preorder numbers subtree(v) as one block [tin[v], tin[v] +
// sz[v]). storing the color at each preorder slot lets add/remove of a subtree
// stream a cache-friendly array -- and it keeps the whole thing iterative, so a
// depth-1e5 path never touches the call stack.
class LomsatGelral {
public:
    // color is 0-indexed by vertex-1 (color[i] = color of vertex i+1); edges are
    // 1-based undirected pairs. returns answers 0-indexed the same way.
    std::vector<std::int64_t> solve(int n, const std::vector<int>& color,
                                    const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        answer_.assign(static_cast<std::size_t>(n) + 1, 0);
        if (n <= 0) return {};

        col_.assign(static_cast<std::size_t>(n) + 1, 0);
        max_color_ = 0;
        for (int v = 1; v <= n; ++v) {
            int c = color[static_cast<std::size_t>(v - 1)];
            col_[static_cast<std::size_t>(v)] = c;
            if (c > max_color_) max_color_ = c;
        }

        build_adjacency(edges);
        flatten();          // iterative preorder -> parent, sz, tin, ord, heavy.
        sack();             // iterative small-to-large over the euler blocks.

        std::vector<std::int64_t> out(static_cast<std::size_t>(n));
        for (int v = 1; v <= n; ++v)
            out[static_cast<std::size_t>(v - 1)] = answer_[static_cast<std::size_t>(v)];
        return out;
    }

private:
    int n_ = 0;
    int max_color_ = 0;                     // largest color, sizes the histogram.
    std::vector<int> col_;                  // col_[v] = color of vertex v, 1-indexed.
    std::vector<std::int64_t> answer_;      // answer_[v], 1-indexed.

    // adjacency as flat CSR: adj_[head_[v] .. head_[v+1]) are v's neighbors.
    // one contiguous block beats n little vectors -- no per-vertex allocation and
    // the walk stays in cache.
    std::vector<int> head_;
    std::vector<int> adj_;

    std::vector<int> parent_;
    std::vector<int> sz_;
    std::vector<int> heavy_;      // heavy_[v] = heaviest child, or 0 if none.
    std::vector<int> tin_;        // preorder index of v, in [0, n).
    std::vector<int> ord_;        // ord_[t] = the vertex at preorder slot t.
    std::vector<int> euler_col_;  // euler_col_[t] = col_[ord_[t]] -- prefolded.

    // the sliding color histogram and its running maximum.
    std::vector<int> cnt_;        // cnt_[color] over the vertices now in the bag.
    int max_cnt_ = 0;             // largest count present.
    std::int64_t sum_at_max_ = 0; // sum of colors that hit max_cnt_.

    // undirected edges -> CSR. count degrees, prefix-sum into head_, scatter.
    void build_adjacency(const std::vector<std::pair<int, int>>& edges) {
        head_.assign(static_cast<std::size_t>(n_) + 2, 0);
        for (const auto& e : edges) {
            ++head_[static_cast<std::size_t>(e.first) + 1];
            ++head_[static_cast<std::size_t>(e.second) + 1];
        }
        for (int v = 1; v <= n_ + 1; ++v)
            head_[static_cast<std::size_t>(v)] += head_[static_cast<std::size_t>(v - 1)];
        adj_.assign(head_[static_cast<std::size_t>(n_ + 1)], 0);
        std::vector<int> fill = head_;  // cursor per vertex, consumed as we scatter.
        for (const auto& e : edges) {
            adj_[static_cast<std::size_t>(fill[static_cast<std::size_t>(e.first)]++)] = e.second;
            adj_[static_cast<std::size_t>(fill[static_cast<std::size_t>(e.second)]++)] = e.first;
        }
    }

    // iterative preorder from root 1. records parent and the preorder order; the
    // reverse of that order is a valid post-order for the size roll-up.
    void flatten() {
        parent_.assign(static_cast<std::size_t>(n_) + 1, 0);
        sz_.assign(static_cast<std::size_t>(n_) + 1, 1);
        sz_[0] = 0;  // vertex 0 is the "no child yet" sentinel -- size zero so the
                     // first real child, even a leaf, wins the heavy comparison.
        heavy_.assign(static_cast<std::size_t>(n_) + 1, 0);
        tin_.assign(static_cast<std::size_t>(n_) + 1, 0);
        ord_.assign(static_cast<std::size_t>(n_), 0);
        euler_col_.assign(static_cast<std::size_t>(n_), 0);

        std::vector<int> stack;
        stack.reserve(static_cast<std::size_t>(n_));
        stack.push_back(1);
        parent_[1] = 0;
        int timer = 0;
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            tin_[static_cast<std::size_t>(v)] = timer;
            ord_[static_cast<std::size_t>(timer)] = v;
            euler_col_[static_cast<std::size_t>(timer)] = col_[static_cast<std::size_t>(v)];
            ++timer;
            for (int i = head_[static_cast<std::size_t>(v)];
                 i < head_[static_cast<std::size_t>(v + 1)]; ++i) {
                int u = adj_[static_cast<std::size_t>(i)];
                if (u == parent_[static_cast<std::size_t>(v)]) continue;
                parent_[static_cast<std::size_t>(u)] = v;
                stack.push_back(u);
            }
        }

        // reverse preorder: a child always precedes its parent, so sizes roll up
        // in one pass and each parent learns its heaviest child on the way.
        for (int t = n_ - 1; t >= 1; --t) {
            int v = ord_[static_cast<std::size_t>(t)];
            int p = parent_[static_cast<std::size_t>(v)];
            sz_[static_cast<std::size_t>(p)] += sz_[static_cast<std::size_t>(v)];
            if (sz_[static_cast<std::size_t>(v)] >
                sz_[static_cast<std::size_t>(heavy_[static_cast<std::size_t>(p)])])
                heavy_[static_cast<std::size_t>(p)] = v;
        }
    }

    // one color enters the bag. the max only ever climbs by one per insert, so a
    // single compare keeps max_cnt_ / sum_at_max_ exact.
    inline void add(int color) {
        int c = ++cnt_[static_cast<std::size_t>(color)];
        if (c > max_cnt_) {
            max_cnt_ = c;
            sum_at_max_ = color;
        } else if (c == max_cnt_) {
            sum_at_max_ += color;
        }
    }

    // fold every vertex of subtree(v) into the bag -- one contiguous euler scan.
    inline void add_subtree(int v) {
        int lo = tin_[static_cast<std::size_t>(v)];
        int hi = lo + sz_[static_cast<std::size_t>(v)];
        for (int t = lo; t < hi; ++t) add(euler_col_[static_cast<std::size_t>(t)]);
    }

    // drop subtree(v) from the bag. this is only ever called when the bag holds
    // EXACTLY subtree(v) -- a light vertex clearing itself -- so afterward every
    // count is zero and the max resets flat. no per-remove max bookkeeping needed,
    // which is the whole reason the aggregate stays a single compare on add.
    inline void remove_subtree(int v) {
        int lo = tin_[static_cast<std::size_t>(v)];
        int hi = lo + sz_[static_cast<std::size_t>(v)];
        for (int t = lo; t < hi; ++t)
            --cnt_[static_cast<std::size_t>(euler_col_[static_cast<std::size_t>(t)])];
        max_cnt_ = 0;
        sum_at_max_ = 0;
    }

    void sack() {
        // histogram is indexed by color, not by vertex -- the problem promises
        // c_v <= n, but sizing to the actual max keeps it correct even if a caller
        // hands over a color past n.
        cnt_.assign(static_cast<std::size_t>(max_color_) + 1, 0);
        max_cnt_ = 0;
        sum_at_max_ = 0;
        heavy_pushed_.assign(static_cast<std::size_t>(n_) + 1, 0);

        // explicit dfs. child_ptr_[v] walks v's neighbor block; we skip the parent
        // and defer the heavy child to last so its bag survives into process(v).
        std::vector<int> child_ptr(head_);  // per-vertex cursor into adj_.
        std::vector<int> stack;
        stack.reserve(static_cast<std::size_t>(n_));
        stack.push_back(1);

        while (!stack.empty()) {
            int v = stack.back();
            int& ci = child_ptr[static_cast<std::size_t>(v)];
            int end = head_[static_cast<std::size_t>(v + 1)];

            // advance to the next light child, if any.
            int next = 0;
            while (ci < end) {
                int u = adj_[static_cast<std::size_t>(ci++)];
                if (u == parent_[static_cast<std::size_t>(v)]) continue;
                if (u == heavy_[static_cast<std::size_t>(v)]) continue;  // heavy last.
                next = u;
                break;
            }
            if (next != 0) {
                stack.push_back(next);  // descend, it will clear itself on return.
                continue;
            }
            // light children exhausted. push the heavy child exactly once, and
            // only after it we fall through to the merge.
            if (ci == end && heavy_[static_cast<std::size_t>(v)] != 0 &&
                !heavy_done(v)) {
                mark_heavy_done(v);
                stack.push_back(heavy_[static_cast<std::size_t>(v)]);
                continue;
            }

            // MERGE. cnt_ already holds the heavy subtree (kept live). add v and
            // each light subtree to complete subtree(v), read the answer, then if
            // v is light wipe the bag for its next sibling.
            add(col_[static_cast<std::size_t>(v)]);
            for (int i = head_[static_cast<std::size_t>(v)]; i < end; ++i) {
                int u = adj_[static_cast<std::size_t>(i)];
                if (u == parent_[static_cast<std::size_t>(v)]) continue;
                if (u == heavy_[static_cast<std::size_t>(v)]) continue;
                add_subtree(u);
            }
            answer_[static_cast<std::size_t>(v)] = sum_at_max_;

            bool keep = (v == 1) ||
                        (v == heavy_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(v)])]);
            if (!keep) remove_subtree(v);
            stack.pop_back();
        }
    }

    // the heavy child must be visited once, and only after the light children.
    // a per-vertex flag records that we already pushed it, so re-inspecting the
    // frame after it returns falls through to the merge instead of pushing it
    // again.
    std::vector<char> heavy_pushed_;
    bool heavy_done(int v) const {
        return heavy_pushed_[static_cast<std::size_t>(v)] != 0;
    }
    void mark_heavy_done(int v) { heavy_pushed_[static_cast<std::size_t>(v)] = 1; }
};

// answers for vertices 1..n. color 0-indexed by vertex-1, edges 1-based undirected,
// root is vertex 1. O(n log n) time, O(n) space.
inline std::vector<std::int64_t> solve(
    int n, const std::vector<int>& color,
    const std::vector<std::pair<int, int>>& edges) {
    LomsatGelral g;
    return g.solve(n, color, edges);
}

}  // namespace p0013
