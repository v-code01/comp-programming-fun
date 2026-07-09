#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace p0039 {

// offline dynamic connectivity. all q ops are known up front -- that's the whole
// game. online you'd reach for a link-cut tree or Euler-tour forest; offline you
// don't have to. you can look at the future.
//
// THE IDEA -- an edge lives on time INTERVALS. walk the op stream once. every add
// opens an interval, the matching remove closes it, and an edge still present at
// the end runs to time q. so each stretch of an edge's life is a half-open
// interval [add, remove) on the leaf axis [0, q). an edge added and removed five
// times is five disjoint intervals -- the fragments are the point.
//
// now hang those intervals on a SEGMENT TREE over the timeline. an interval
// [l, r) decomposes into O(log q) canonical segtree nodes that exactly tile it --
// the same decomposition a range-update segtree uses. push the edge onto each of
// those nodes. every edge costs O(log q) placements, O(q log q) total.
//
// then DFS the tree. at a node, UNITE all edges parked there in the DSU; those
// edges are alive across the node's whole time span, so they belong in the
// structure for every leaf below it. recurse into both children. at a LEAF -- a
// single instant -- the DSU holds exactly the edges whose intervals cover that
// instant, so a query reads off connected / not. on the way back UP, ROLL BACK the
// unions this node made, so the sibling subtree sees a clean structure. the union
// state is a stack that grows going down and unwinds coming up.
//
// WHY NO PATH COMPRESSION. rollback needs every mutation reversible. path
// compression rewrites parent pointers all along a find -- an unbounded, unlogged
// cascade you can't cheaply undo. so it's out. union by rank alone, no shortcuts.
// and that's enough: uniting the shorter tree under the taller one keeps height
// <= log2(size). a root of rank r subtends at least 2^r nodes -- induction: rank
// only rises when two rank-r roots merge, and that doubles the count. so rank <=
// log2 n, find walks at most log2 n parents, O(log n) worst case with zero
// amortization owed. exactly the property compression would have stolen.
//
// COMPLEXITY. reading the ops is Omega(q). each edge lands on O(log q) nodes and
// each union / find is O(log n) with rank only -- but wrapped in the DSU's
// near-constant behavior it's O(alpha) amortized in practice, O(log n) worst per
// call. total O((n + q log q) * alpha). the DFS touches every placement once going
// down and rolls it back once coming up: linear in placements. no online structure
// beats this offline, and the naive recompute-per-query oracle is O(q * (n + q)) --
// this trades that q for a log q.
//
// DEPTH. the segtree is O(log q) tall -- ~18 levels at q = 2e5. the DFS recursion
// rides that height, not q, so the call stack never gets near trouble. no need for
// an explicit stack; the recursion is the clean way to say it.

// rollback DSU -- union by rank, no path compression, an undo log. every union
// pushes one record; a checkpoint is just the log's length, and rollback pops back
// to it. that's the entire trick that lets the segtree DFS share one structure.
class RollbackDsu {
public:
    // vertices 0..n-1. callers here pass n+1 and use 1..n; slot 0 rides along
    // unused. parent points to self, rank starts flat.
    void init(int n) {
        parent_.resize(static_cast<std::size_t>(n));
        rank_.assign(static_cast<std::size_t>(n), 0);
        for (int i = 0; i < n; ++i) parent_[static_cast<std::size_t>(i)] = i;
        log_.clear();
    }

    // no compression -- walk to the root and leave the path as it was, or rollback
    // couldn't put it back. rank keeps this at O(log n).
    int find(int x) const {
        while (parent_[static_cast<std::size_t>(x)] != x)
            x = parent_[static_cast<std::size_t>(x)];
        return x;
    }

    // hang the shorter tree under the taller. returns false when they already share
    // a root -- nothing changed, nothing logged, so rollback stays exact.
    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        // a is the taller-or-equal root; b sinks under it.
        if (rank_[static_cast<std::size_t>(a)] < rank_[static_cast<std::size_t>(b)])
            std::swap(a, b);
        parent_[static_cast<std::size_t>(b)] = a;
        // rank rises only on a tie -- that's the one case the root got taller.
        bool bumped =
            rank_[static_cast<std::size_t>(a)] == rank_[static_cast<std::size_t>(b)];
        if (bumped) ++rank_[static_cast<std::size_t>(a)];
        log_.push_back({b, bumped});
        return true;
    }

    bool connected(int a, int b) const { return find(a) == find(b); }

    // a checkpoint is the log length -- cheap to take, cheap to compare.
    std::size_t checkpoint() const { return log_.size(); }

    // undo every union back to a checkpoint, newest first. the child still points
    // at the root it was merged into, so read the root from there, cut the link,
    // and un-bump the rank if this union raised it.
    void rollback(std::size_t mark) {
        while (log_.size() > mark) {
            Record r = log_.back();
            log_.pop_back();
            int root = parent_[static_cast<std::size_t>(r.child)];
            parent_[static_cast<std::size_t>(r.child)] = r.child;
            if (r.bumped) --rank_[static_cast<std::size_t>(root)];
        }
    }

private:
    struct Record {
        int child;    // the root that got reparented -- set parent[child]=child to undo.
        bool bumped;  // did this union raise the surviving root's rank?
    };
    std::vector<int> parent_;
    std::vector<int> rank_;      // height bound, not size -- that's what caps find.
    std::vector<Record> log_;    // the undo stack. grows down the DFS, unwinds up.
};

// one op in the statement's alphabet, vertices 1-based.
//   type 1: add edge u--v
//   type 2: remove edge u--v (present by contract)
//   type 3: query connected(u, v) -> 1 / 0
struct Op {
    int type;
    int u;
    int v;
};

// the whole offline pipeline: scan ops into edge-intervals, tile them on a segtree
// over time, DFS with the rollback DSU, and read off the type-3 answers in order.
class OfflineDynamicConnectivity {
public:
    // n vertices, ops in timeline order. returns one 0/1 per type-3 query, in the
    // order the queries appear.
    std::vector<int> solve(int n, const std::vector<Op>& ops) {
        q_ = static_cast<int>(ops.size());
        answers_.clear();
        if (q_ == 0) return answers_;

        // 1..n are real vertices; init n+1 so 1-based indices land in slot.
        dsu_.init(n + 1);

        edges_.clear();
        seg_.assign(static_cast<std::size_t>(4 * q_), {});
        is_query_.assign(static_cast<std::size_t>(q_), 0);
        qu_.assign(static_cast<std::size_t>(q_), 0);
        qv_.assign(static_cast<std::size_t>(q_), 0);

        // walk the stream once, opening and closing intervals. an edge in flight is
        // keyed by its normalized endpoints and remembers when it was added.
        std::unordered_map<std::int64_t, int> open;  // edge key -> add time.
        open.reserve(static_cast<std::size_t>(q_) * 2);

        for (int i = 0; i < q_; ++i) {
            const Op& o = ops[static_cast<std::size_t>(i)];
            if (o.type == 1) {
                open[key(n, o.u, o.v)] = i;  // interval opens at this instant.
            } else if (o.type == 2) {
                auto it = open.find(key(n, o.u, o.v));
                // present by contract -- guarded anyway so a bad input can't corrupt.
                if (it != open.end()) {
                    addInterval(o.u, o.v, it->second, i);
                    open.erase(it);
                }
            } else {
                is_query_[static_cast<std::size_t>(i)] = 1;
                qu_[static_cast<std::size_t>(i)] = o.u;
                qv_[static_cast<std::size_t>(i)] = o.v;
            }
        }
        // edges never removed live to the end of the timeline.
        for (const auto& kv : open) {
            int add_time = kv.second;
            int u = static_cast<int>(kv.first / (n + 1));
            int v = static_cast<int>(kv.first % (n + 1));
            addInterval(u, v, add_time, q_);
        }

        // leaves are visited left to right, so answers land in query order.
        answers_.reserve(static_cast<std::size_t>(q_));
        dfs(1, 0, q_);
        return answers_;
    }

private:
    int q_ = 0;
    RollbackDsu dsu_;
    std::vector<std::pair<int, int>> edges_;    // edge id -> (u, v).
    std::vector<std::vector<int>> seg_;         // segtree node -> edge ids parked there.
    std::vector<char> is_query_;                // leaf -> is it a type-3 op?
    std::vector<int> qu_, qv_;                  // leaf -> query endpoints.
    std::vector<int> answers_;

    // pack an unordered edge into one int64. n <= 2e5 so (n+1) fits the radix and
    // the product stays well inside int64. min goes high, max goes low.
    static std::int64_t key(int n, int u, int v) {
        int a = u < v ? u : v;
        int b = u < v ? v : u;
        return static_cast<std::int64_t>(a) * (n + 1) + b;
    }

    // record an edge alive over [l, r), then tile that interval onto the segtree.
    void addInterval(int u, int v, int l, int r) {
        if (l >= r) return;  // zero-width -- add and remove at the same instant.
        int e = static_cast<int>(edges_.size());
        edges_.push_back({u, v});
        add(1, 0, q_, l, r, e);
    }

    // canonical range decomposition -- push e onto the O(log q) nodes that tile
    // [l, r). same recursion a lazy-propagation segtree uses, minus the laziness.
    void add(int node, int lo, int hi, int l, int r, int e) {
        if (r <= lo || hi <= l) return;
        if (l <= lo && hi <= r) {
            seg_[static_cast<std::size_t>(node)].push_back(e);
            return;
        }
        int mid = (lo + hi) / 2;
        add(node * 2, lo, mid, l, r, e);
        add(node * 2 + 1, mid, hi, l, r, e);
    }

    // the descent. unite this node's edges, checkpoint before so the exact set can
    // be undone, recurse, answer at leaves, roll back on the way out. depth is the
    // tree height ~ log q, never q.
    void dfs(int node, int lo, int hi) {
        std::size_t mark = dsu_.checkpoint();
        for (int e : seg_[static_cast<std::size_t>(node)])
            dsu_.unite(edges_[static_cast<std::size_t>(e)].first,
                       edges_[static_cast<std::size_t>(e)].second);

        if (hi - lo == 1) {
            // a single instant. the DSU is exactly the graph at time lo.
            if (is_query_[static_cast<std::size_t>(lo)]) {
                answers_.push_back(dsu_.connected(qu_[static_cast<std::size_t>(lo)],
                                                  qv_[static_cast<std::size_t>(lo)])
                                       ? 1
                                       : 0);
            }
        } else {
            int mid = (lo + hi) / 2;
            dfs(node * 2, lo, mid);
            dfs(node * 2 + 1, mid, hi);
        }

        dsu_.rollback(mark);  // leave the structure as the sibling found it.
    }
};

// convenience wrapper -- one call, one answer stream.
inline std::vector<int> solve(int n, const std::vector<Op>& ops) {
    OfflineDynamicConnectivity solver;
    return solver.solve(n, ops);
}

}  // namespace p0039
