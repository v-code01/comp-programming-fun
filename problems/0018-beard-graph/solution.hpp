#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0018 {

// Beard Graph -- a tree whose n-1 edges each carry a color. every edge starts
// black. flip an edge white, flip it back black, or ask about the simple path
// between two vertices: if every edge on it is black, the answer is the number
// of edges (the path's edge-distance); if any edge is white, the answer is -1;
// if the endpoints coincide, 0.
//
// the naive read walks the a->b path per query -- O(n) each, O(nm) at the
// limits, 1e10 operations. dead. the trick is HEAVY-LIGHT DECOMPOSITION: cut the
// tree into O(log n) heavy chains so any root-to-node path crosses at most
// O(log n) of them, lay every chain down contiguously in one array, and answer a
// path query as O(log n) range queries on a segment tree over that array.
//
// EDGE-ON-DEEPER-VERTEX. the segment tree indexes vertices, but the colors live
// on edges -- one fewer of them. root the tree at vertex 0. every non-root vertex
// v has exactly one parent edge, so pin that edge's color to v's position. the
// root owns no edge; its position is a neutral (len=0, white=0) that never
// counts. this bijection (edge <-> its deeper endpoint) is what lets an
// edge-colored tree ride a vertex-indexed segment tree.
//
// SKIP THE LCA. a path a->b splits at L = LCA(a, b) into a->L and b->L. every
// edge on it hangs below some vertex on those two legs -- except one position we
// must NOT touch: L's own position, which holds the edge ABOVE L, and that edge
// is not on the path. so the same-chain final step aggregates the half-open
// (pos[L], pos[deeper]] -- positions pos[L]+1 .. pos[deeper], L excluded. every
// chain jump above L aggregates [pos[head], pos[node]] whole, because the edge
// above each chain head IS on the path (we hop to parent[head] next). getting
// this one exclusion wrong is the classic HLD-on-edges bug; the differential vs
// the path-walk oracle is what pins it down.
//
// WHAT A POSITION STORES. for the edge at vertex v: len = 1 if black else 0,
// white = 1 if white else 0. a path query sums both across its O(log n) ranges.
// white > 0 means some edge on the path is white -> -1. otherwise every edge is
// black, so the black count len IS the edge count -> that's the answer. a paint
// is one point update flipping a single position between (1,0) and (0,1).
//
// COMPLEXITY. build is O(n): two linear sweeps (sizes + heavy child, then chain
// assignment) plus an O(n) segment-tree build. a paint is O(log n) -- one point
// update. a type-3 query is O(log^2 n): O(log n) chains, each an O(log n)
// segment-tree range query. whole thing O((n + m) log^2 n).
//
// LOWER BOUND, honest. you cannot answer without reading the n-1 edges and the m
// queries, so Omega(n + m) is forced -- every edge and every query has to be
// looked at once. HLD pays a log^2 n multiplier on the query mix: O(log n) chains
// times O(log n) per segment-tree range. no claim this is optimal in the absolute
// -- link/cut trees do the same job in O(log n) amortized per op, one log better
// -- only that O((n + m) log^2 n) sits two logs off the read-the-input floor, and
// both logs are structural: chain count, then range-query depth.
//
// EVERYTHING iterative. n reaches 1e5 and the tree can be a single path, so a
// recursive size DFS or a recursive chain build would blow the stack. the two
// build sweeps run off an explicit preorder and a reverse of it; the segment tree
// is a bottom-up array. a path of 1e5 vertices is a long loop, never a deep frame
// stack.
class BeardHLD {
public:
    // build over an n-vertex tree from 0-based edges, edge i = edges[i]. every
    // edge starts black -- the statement's initial condition, folded in so the
    // caller never re-asserts it. n-1 edges, connected: the caller's contract.
    void build(int n, const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        adj_.assign(static_cast<std::size_t>(n), {});
        // remember, per edge, which endpoints it joins -- after rooting, the
        // deeper one owns the edge's segment position.
        edge_ends_ = edges;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto [u, v] = edges[i];
            adj_[static_cast<std::size_t>(u)].push_back(v);
            adj_[static_cast<std::size_t>(v)].push_back(u);
        }

        parent_.assign(static_cast<std::size_t>(n), -1);
        depth_.assign(static_cast<std::size_t>(n), 0);
        heavy_.assign(static_cast<std::size_t>(n), -1);
        size_.assign(static_cast<std::size_t>(n), 0);
        head_.assign(static_cast<std::size_t>(n), 0);
        pos_.assign(static_cast<std::size_t>(n), 0);
        edge_pos_.assign(edges.size(), 0);

        root_pass();   // parent, depth, subtree sizes, heavy child.
        chain_pass();  // head[], pos[], contiguous heavy chains.

        // map each input edge to the position of its deeper endpoint, and seed the
        // segment tree: every real edge black (1, 0), the root's slot neutral.
        seg_len_.assign(2 * static_cast<std::size_t>(n), 0);
        seg_white_.assign(2 * static_cast<std::size_t>(n), 0);
        for (std::size_t i = 0; i < edges.size(); ++i) {
            const auto [u, v] = edge_ends_[i];
            const int child = (parent_[static_cast<std::size_t>(u)] == v) ? u : v;
            const int p = pos_[static_cast<std::size_t>(child)];
            edge_pos_[i] = p;
            seg_len_[static_cast<std::size_t>(n_ + p)] = 1;  // black.
        }
        for (int i = n_ - 1; i >= 1; --i) pull(i);
    }

    // paint edge i black -- its position becomes (len=1, white=0). idempotent.
    void paint_black(int i) { set_edge(i, 1, 0); }

    // paint edge i white -- its position becomes (len=0, white=1). idempotent.
    void paint_white(int i) { set_edge(i, 0, 1); }

    // the simple path a -> b. all edges black -> edge count; any white -> -1;
    // a == b -> 0. a, b are 0-based vertices.
    int query(int a, int b) const {
        if (a == b) return 0;  // no edges on a zero-length path.
        std::int64_t total_len = 0;
        std::int64_t total_white = 0;

        // climb the two endpoints chain by chain until they share a chain. the
        // deeper-head endpoint jumps; its chain segment [pos[head], pos[u]] is
        // fully on the path, edge above head included.
        while (head_[static_cast<std::size_t>(a)] !=
               head_[static_cast<std::size_t>(b)]) {
            if (depth_[static_cast<std::size_t>(head_[static_cast<std::size_t>(a)])] <
                depth_[static_cast<std::size_t>(head_[static_cast<std::size_t>(b)])]) {
                std::swap(a, b);
            }
            const int h = head_[static_cast<std::size_t>(a)];
            accumulate(pos_[static_cast<std::size_t>(h)],
                       pos_[static_cast<std::size_t>(a)], total_len, total_white);
            a = parent_[static_cast<std::size_t>(h)];
        }

        // same chain now. the shallower vertex is L = LCA(a, b). its own position
        // holds the edge ABOVE L -- not on the path -- so aggregate strictly below
        // it: pos[L]+1 .. pos[deeper]. skipping pos[L] is the whole edge trick.
        if (depth_[static_cast<std::size_t>(a)] >
            depth_[static_cast<std::size_t>(b)]) {
            std::swap(a, b);
        }
        if (pos_[static_cast<std::size_t>(a)] + 1 <=
            pos_[static_cast<std::size_t>(b)]) {
            accumulate(pos_[static_cast<std::size_t>(a)] + 1,
                       pos_[static_cast<std::size_t>(b)], total_len, total_white);
        }

        if (total_white > 0) return -1;
        return static_cast<int>(total_len);  // all black -> len is the edge count.
    }

private:
    // ---- build sweep 1: root the tree, size it, pick heavy children ----------
    // one explicit preorder off a stack (no recursion), then subtree sizes in the
    // reverse of that order -- a child always precedes its parent in the reverse,
    // so its size is final before the parent reads it.
    void root_pass() {
        order_.clear();
        order_.reserve(static_cast<std::size_t>(n_));
        std::vector<int> stack;
        stack.reserve(static_cast<std::size_t>(n_));
        stack.push_back(0);
        parent_[0] = -1;
        depth_[0] = 0;
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            order_.push_back(v);
            for (const int w : adj_[static_cast<std::size_t>(v)]) {
                if (w == parent_[static_cast<std::size_t>(v)]) continue;
                parent_[static_cast<std::size_t>(w)] = v;
                depth_[static_cast<std::size_t>(w)] =
                    depth_[static_cast<std::size_t>(v)] + 1;
                stack.push_back(w);
            }
        }
        for (std::size_t i = order_.size(); i-- > 0;) {
            const int v = order_[i];
            size_[static_cast<std::size_t>(v)] += 1;  // count self.
            int best = 0;
            for (const int w : adj_[static_cast<std::size_t>(v)]) {
                if (w == parent_[static_cast<std::size_t>(v)]) continue;
                size_[static_cast<std::size_t>(v)] += size_[static_cast<std::size_t>(w)];
                if (size_[static_cast<std::size_t>(w)] > best) {
                    best = size_[static_cast<std::size_t>(w)];
                    heavy_[static_cast<std::size_t>(v)] = w;  // heaviest child.
                }
            }
        }
    }

    // ---- build sweep 2: lay chains down contiguously -------------------------
    // pop a chain head, walk its heavy spine assigning consecutive positions, and
    // push each light child as the head of its own future chain. because the whole
    // spine is numbered before any light subtree, every chain occupies one
    // contiguous position range with the head smallest -- exactly what the range
    // query [pos[head], pos[node]] needs. no recursion: the stack holds pending
    // chain heads, never a frame per tree level.
    void chain_pass() {
        int timer = 0;
        std::vector<std::pair<int, int>> stack;  // (chain start, its head).
        stack.reserve(static_cast<std::size_t>(n_));
        stack.emplace_back(0, 0);
        while (!stack.empty()) {
            auto [v, h] = stack.back();
            stack.pop_back();
            while (v != -1) {
                head_[static_cast<std::size_t>(v)] = h;
                pos_[static_cast<std::size_t>(v)] = timer++;
                for (const int w : adj_[static_cast<std::size_t>(v)]) {
                    if (w == parent_[static_cast<std::size_t>(v)]) continue;
                    if (w == heavy_[static_cast<std::size_t>(v)]) continue;
                    stack.emplace_back(w, w);  // light child starts a new chain.
                }
                v = heavy_[static_cast<std::size_t>(v)];  // stay on the spine.
            }
        }
    }

    // ---- segment tree: bottom-up, point update + range sum of (len, white) ---
    void pull(int node) {
        seg_len_[static_cast<std::size_t>(node)] =
            seg_len_[static_cast<std::size_t>(2 * node)] +
            seg_len_[static_cast<std::size_t>(2 * node + 1)];
        seg_white_[static_cast<std::size_t>(node)] =
            seg_white_[static_cast<std::size_t>(2 * node)] +
            seg_white_[static_cast<std::size_t>(2 * node + 1)];
    }

    // overwrite edge i's leaf with (lv, wv) and mend the path to the root.
    void set_edge(int i, int lv, int wv) {
        int node = n_ + edge_pos_[static_cast<std::size_t>(i)];
        seg_len_[static_cast<std::size_t>(node)] = lv;
        seg_white_[static_cast<std::size_t>(node)] = wv;
        for (node >>= 1; node >= 1; node >>= 1) pull(node);
    }

    // sum (len, white) over the inclusive position range [l, r] into the running
    // totals. standard half-open bottom-up walk on [l, r+1).
    void accumulate(int l, int r, std::int64_t& tl, std::int64_t& tw) const {
        for (l += n_, r += n_ + 1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) {
                tl += seg_len_[static_cast<std::size_t>(l)];
                tw += seg_white_[static_cast<std::size_t>(l)];
                ++l;
            }
            if (r & 1) {
                --r;
                tl += seg_len_[static_cast<std::size_t>(r)];
                tw += seg_white_[static_cast<std::size_t>(r)];
            }
        }
    }

    int n_ = 0;
    std::vector<std::vector<int>> adj_;
    std::vector<std::pair<int, int>> edge_ends_;  // endpoints of each input edge.
    std::vector<int> parent_;
    std::vector<int> depth_;
    std::vector<int> heavy_;  // heaviest child, or -1 for a leaf.
    std::vector<int> size_;   // subtree size, rooted at 0.
    std::vector<int> head_;   // chain head of each vertex.
    std::vector<int> pos_;    // segment-tree position of each vertex.
    std::vector<int> order_;  // preorder, reused for the reverse size sweep.
    std::vector<int> edge_pos_;      // input edge i -> its deeper endpoint's pos.
    std::vector<int> seg_len_;       // segment tree: black-edge count per node.
    std::vector<int> seg_white_;     // segment tree: white-edge count per node.
};

}  // namespace p0018
