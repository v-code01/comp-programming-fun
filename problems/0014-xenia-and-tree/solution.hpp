#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace p0014 {

// Xenia and Tree -- nearest red vertex on an unweighted tree, online.
// vertex 1 (index 0 here) starts red. paint a vertex red, or ask how many edges
// separate a vertex from the closest red one. answers are exact, no waiting for
// the query stream to end.
//
// the naive read of the problem is a BFS per query -- O(n) each, O(nm) total,
// 1e10 for the limits. dead. the trick is CENTROID DECOMPOSITION: a second tree
// on the same vertices, O(log n) tall, that lets every paint and every query
// touch just O(log n) vertices.
//
// THE CENTROID TREE. pick a centroid of the whole tree -- a vertex whose removal
// leaves every piece with at most half the vertices. make it the root. remove
// it, recurse on each piece, and hang each piece's centroid under it as a child.
// the halving is the whole point: a vertex sits in at most O(log n) of these
// nested pieces, so it has at most O(log n) centroid ancestors. that bound is the
// budget every operation spends.
//
// THE INVARIANT best[c]. for each centroid c we keep
//     best[c] = min over red vertices r inside c's centroid-subtree of dist(c, r)
// where dist is the ordinary edge count in the ORIGINAL tree. start it at +inf.
// c's centroid-subtree is exactly the piece c was the centroid of -- the set of
// vertices that had c as an ancestor when the decomposition split them off.
//
// PAINT v. v belongs to the centroid-subtree of every one of its centroid
// ancestors c (that's what "ancestor" means here). so painting v red can only
// improve best[c] for those c, and no others:
//     for each centroid ancestor c of v:  best[c] = min(best[c], dist(c, v))
// O(log n) ancestors, one relaxation each.
//
// QUERY v. claim: the answer is
//     min over centroid ancestors c of v of ( dist(c, v) + best[c] ).
// WHY it's an upper bound -- every term is a real path. best[c] = dist(c, r) for
// some red r, and dist(c, v) + dist(c, r) is the length of a v->c->r walk, which
// is >= dist(v, r) >= the true answer. so no term ever undershoots the truth.
//
// WHY it's tight -- here's the load-bearing fact about centroid trees. let r* be
// the genuinely nearest red vertex and let c* = LCA(v, r*) in the CENTROID tree.
// c* is a centroid ancestor of both v and r*. and c* lies on the original-tree
// path from v to r*: at the level where c* was chosen, v and r* were in the same
// piece, and c* is the centroid whose removal split that piece -- v and r* land
// in different fragments (or one of them IS c*), so every original-tree path
// between them must cross c*. crossing c* means
//     dist(v, r*) = dist(v, c*) + dist(c*, r*).
// r* is red and inside c*'s subtree, so best[c*] <= dist(c*, r*), giving
//     dist(v, c*) + best[c*] <= dist(v, c*) + dist(c*, r*) = dist(v, r*).
// the min over ancestors is <= this <= the true answer, and by the upper bound it
// can't be less. equal. the min over c of (dist(c,v) + best[c]) is the answer,
// every time.
//
// DISTANCES, precomputed. building the decomposition, the moment we fix a
// centroid c we BFS its piece and stamp dist(c, u) onto every u in it -- that's
// exactly u gaining c as an ancestor. we store the pair (c, dist(c,u)) on u. a
// vertex collects one pair per ancestor, O(log n) of them, O(n log n) pairs total
// (~1.7M at n=1e5). paint and query then just walk u's own short list -- no LCA
// machinery, no per-op tree walk, the ancestors are sitting right there.
//
// COMPLEXITY. build is O(n log n): each of the O(log n) levels does O(n) total
// work across its pieces (size counts, centroid descents, distance BFS). each
// paint and each query is O(log n) -- the length of one vertex's ancestor list.
// whole thing O((n + m) log n).
//
// LOWER BOUND, honest. you cannot answer without reading the tree and the query
// stream, so Omega(n + m) is forced -- n-1 edges and m operations, each has to be
// looked at. centroid decomposition pays a log n multiplier on top: every update
// and every query touches the O(log n) centroid ancestors, and that depth is
// intrinsic to this technique. no claim it's optimal in the absolute -- link/cut
// or Euler-tour tricks exist -- only that O((n + m) log n) is a log factor off the
// read-the-input floor, and the log is the centroid-tree height.
//
// EVERYTHING iterative. n reaches 1e5 and the tree can be a path, so a recursive
// size count or a recursive centroid build would blow the stack. size counting,
// centroid descent, and the distance BFS are all explicit-queue loops -- a path
// of 1e5 vertices is just a long queue, never a deep frame stack.
class NearestRed {
public:
    // build the centroid tree over an n-vertex tree from 0-based edges.
    // n-1 edges, connected -- caller's contract, the statement guarantees it.
    // vertex 0 (vertex 1 in the statement) is painted red here, folded in so the
    // caller never has to remember the initial condition.
    void build(int n, const std::vector<std::pair<int, int>>& edges) {
        n_ = n;
        adj_.assign(static_cast<std::size_t>(n), {});
        for (const auto& [u, v] : edges) {
            adj_[static_cast<std::size_t>(u)].push_back(v);
            adj_[static_cast<std::size_t>(v)].push_back(u);
        }

        removed_.assign(static_cast<std::size_t>(n), 0);
        sub_.assign(static_cast<std::size_t>(n), 0);
        par_.assign(static_cast<std::size_t>(n), -1);
        seen_.assign(static_cast<std::size_t>(n), 0);
        dist_.assign(static_cast<std::size_t>(n), 0);
        anc_.assign(static_cast<std::size_t>(n), {});
        best_.assign(static_cast<std::size_t>(n), kInf);
        stamp_ = 0;

        decompose();

        paint(0);  // the statement's initial red vertex.
    }

    // paint v red. relax best[] on each of v's O(log n) centroid ancestors.
    // idempotent -- repainting a red vertex relaxes against distances already
    // recorded and changes nothing.
    void paint(int v) {
        for (const auto& [c, d] : anc_[static_cast<std::size_t>(v)]) {
            if (d < best_[static_cast<std::size_t>(c)]) {
                best_[static_cast<std::size_t>(c)] = d;
            }
        }
    }

    // edges to the nearest red vertex. vertex 0 is always red, so on a connected
    // tree this always finds something and the sentinel never escapes.
    int query(int v) const {
        int ans = kInf;
        for (const auto& [c, d] : anc_[static_cast<std::size_t>(v)]) {
            const int cand = d + best_[static_cast<std::size_t>(c)];
            if (cand < ans) ans = cand;
        }
        return ans;
    }

private:
    // 1<<29 leaves headroom: the largest real distance is < n <= 1e5, so
    // dist + kInf stays well under INT_MAX and never wraps in the query min.
    static constexpr int kInf = 1 << 29;

    // drive the whole decomposition off one work list of (component seed). no
    // recursion -- a component is represented by any one of its vertices, and we
    // never carry a call frame per level.
    void decompose() {
        std::vector<int> seeds;
        seeds.reserve(static_cast<std::size_t>(n_));
        seeds.push_back(0);

        std::vector<int> order;  // BFS order within a component, reused.
        order.reserve(static_cast<std::size_t>(n_));

        while (!seeds.empty()) {
            const int start = seeds.back();
            seeds.pop_back();
            if (removed_[static_cast<std::size_t>(start)]) continue;  // defensive.

            const int comp_size = collect(start, order);
            const int c = find_centroid(start, comp_size, order);

            record_distances(c);  // stamp dist(c, u) onto every u in the piece.

            removed_[static_cast<std::size_t>(c)] = 1;
            for (const int w : adj_[static_cast<std::size_t>(c)]) {
                if (!removed_[static_cast<std::size_t>(w)]) seeds.push_back(w);
            }
        }
    }

    // BFS the live component containing start: fill `order`, set par_ (rooted at
    // start), then accumulate sub_ sizes in reverse BFS order. returns the size.
    // a fresh stamp marks membership so no O(n) clear is needed between calls.
    int collect(int start, std::vector<int>& order) {
        ++stamp_;
        order.clear();
        seen_[static_cast<std::size_t>(start)] = stamp_;
        par_[static_cast<std::size_t>(start)] = -1;
        order.push_back(start);

        for (std::size_t i = 0; i < order.size(); ++i) {
            const int u = order[i];
            for (const int w : adj_[static_cast<std::size_t>(u)]) {
                if (removed_[static_cast<std::size_t>(w)]) continue;
                if (seen_[static_cast<std::size_t>(w)] == stamp_) continue;
                seen_[static_cast<std::size_t>(w)] = stamp_;
                par_[static_cast<std::size_t>(w)] = u;
                order.push_back(w);
            }
        }

        for (const int u : order) sub_[static_cast<std::size_t>(u)] = 1;
        for (std::size_t i = order.size(); i-- > 0;) {
            const int u = order[i];
            const int p = par_[static_cast<std::size_t>(u)];
            if (p != -1) sub_[static_cast<std::size_t>(p)] += sub_[static_cast<std::size_t>(u)];
        }
        return static_cast<int>(order.size());
    }

    // descend from start toward the heavy side until no neighbor holds more than
    // half the component -- that vertex is the centroid. sub_/par_ come from
    // collect(); we came from `prev` (always the parent side, already known
    // light), so we only weigh the child sides.
    int find_centroid(int start, int comp_size, const std::vector<int>& /*order*/) {
        int c = start;
        int prev = -1;
        while (true) {
            int next = -1;
            for (const int w : adj_[static_cast<std::size_t>(c)]) {
                if (removed_[static_cast<std::size_t>(w)] || w == prev) continue;
                // w's side-size seen from c: a child contributes its subtree; the
                // one parent-side neighbor (== prev during descent) would be the
                // rest, but prev is skipped, so this reduces to the child case.
                const int side =
                    (par_[static_cast<std::size_t>(w)] == c)
                        ? sub_[static_cast<std::size_t>(w)]
                        : comp_size - sub_[static_cast<std::size_t>(c)];
                if (side * 2 > comp_size) {
                    next = w;
                    break;
                }
            }
            if (next == -1) break;
            prev = c;
            c = next;
        }
        return c;
    }

    // BFS the piece from its centroid c (not yet removed, so the BFS spans exactly
    // the live component) and append (c, dist(c, u)) to every reached u -- u thus
    // gains c as a centroid ancestor. c itself gets (c, 0).
    void record_distances(int c) {
        ++stamp_;
        std::vector<int>& q = bfs_;
        q.clear();
        seen_[static_cast<std::size_t>(c)] = stamp_;
        dist_[static_cast<std::size_t>(c)] = 0;
        q.push_back(c);
        anc_[static_cast<std::size_t>(c)].emplace_back(c, 0);

        for (std::size_t i = 0; i < q.size(); ++i) {
            const int u = q[i];
            const int du = dist_[static_cast<std::size_t>(u)];
            for (const int w : adj_[static_cast<std::size_t>(u)]) {
                if (removed_[static_cast<std::size_t>(w)]) continue;
                if (seen_[static_cast<std::size_t>(w)] == stamp_) continue;
                seen_[static_cast<std::size_t>(w)] = stamp_;
                dist_[static_cast<std::size_t>(w)] = du + 1;
                q.push_back(w);
                anc_[static_cast<std::size_t>(w)].emplace_back(c, du + 1);
            }
        }
    }

    int n_ = 0;
    int stamp_ = 0;
    std::vector<std::vector<int>> adj_;
    std::vector<char> removed_;
    std::vector<int> sub_;   // subtree sizes, rooted per collect() call.
    std::vector<int> par_;   // rooted parent, per collect() call.
    std::vector<int> seen_;  // membership stamp, avoids O(n) clears.
    std::vector<int> dist_;  // dist from current centroid, per record_distances.
    std::vector<int> bfs_;   // distance-BFS queue, reused.
    // anc_[v] = { (centroid ancestor c, dist(c, v)) }, O(log n) entries.
    std::vector<std::vector<std::pair<int, int>>> anc_;
    std::vector<int> best_;  // best[c] = min dist(c, red) in c's subtree, +inf.
};

}  // namespace p0014
