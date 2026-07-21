#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace p0046 {

// fully dynamic connectivity, ONLINE. n vertices, a stream of ops answered in
// order, no lookahead:
//   insert {u,v}   -- add an undirected edge (absent right now).
//   delete {u,v}   -- remove one (present right now).
//   connected(u,v) -- are they in the same component this instant?
//
// offline you'd sort the ops by time, hang each edge's life-interval on a segment
// tree, and sweep with a rollback DSU -- O(log q) per op, dead simple. see 0039.
// but that reads the WHOLE stream first. online forbids it: you must answer query
// k knowing only ops 0..k. that one restriction is the entire reason this file is
// hard. no interval to hang, no future to peek at -- so we carry a real spanning
// forest and repair it on every delete.
//
// THE STRUCTURE -- Holm, de Lichtenberg, Thorup (2001). keep a spanning forest of
// the graph. every edge, tree or not, holds a LEVEL in [0, L], L = ceil(log2 n).
// the levels obey one invariant that pays for everything:
//
//     a level-i tree edge sits in a tree whose vertex count is <= 2^(L - i).
//
// so a level-0 tree can be the whole graph; a level-L tree is a single vertex.
// levels only ever RISE. that monotonicity is the potential: each of the m edges
// climbs at most L times, so the total work spent hunting replacements is
// O(m log^2 n), i.e. O(log^2 n) amortized per op.
//
// EULER-TOUR TREES. each tree in each level's forest is stored as its Euler tour
// in a balanced BST (a treap keyed by position). the tour is a closed walk that
// crosses every tree edge exactly twice -- once per direction -- so an edge {u,v}
// contributes two "arc" tokens, and every vertex contributes one stable "vertex"
// node. connectivity is "same treap root"; link splices two tours; cut carves the
// sub-tour bounded by an edge's two tokens out of the middle. all O(log n). the
// vertex node is stable across rerootings, which is where we hang the aggregate
// bits the replacement search descends by:
//   - avm: does this subtree hold a vertex with an incident NON-tree edge at this
//          level?  (candidates to reconnect with.)
//   - aem: does this subtree hold a tree edge OF this exact level?  (candidates to
//          push up.)
//
// LEVELS, PLURAL. one forest can't answer "the level-i tree containing u" -- that's
// a subgraph, not the whole component. so keep L+1 forests, F_0 superset F_1 ...
// superset F_L, where F_i spans exactly the edges of level >= i. F_0 is the full
// spanning forest and owns every connectivity query. an edge of level l lives as
// tree tokens in F_0..F_l. rising a level = splicing the edge into one more forest.
//
// DELETE, the whole game. deleting a non-tree edge is bookkeeping -- drop it from
// its two incidence lists. deleting a TREE edge e of level l cuts it out of
// F_0..F_l, splitting each into two, and then hunts a replacement, level by level
// from l down to 0:
//   1. of the two split pieces in F_i, take the SMALLER, T_small. it has <=
//      2^(L-i-1) vertices -- half the pre-split tree -- which is exactly 2^(L-(i+1)).
//   2. push every level-i tree edge inside T_small up to level i+1. legal precisely
//      because of the bound in step 1: their new level-(i+1) tree still fits under
//      2^(L-(i+1)). this is where the amortization is banked.
//   3. scan T_small's level-i NON-tree edges. the first one whose far endpoint is
//      in the OTHER piece is the replacement -- promote it to a tree edge at level i
//      and stop, the forest is whole again. any that stays inside T_small can't help
//      at this level, so push it up to i+1 too and keep scanning.
// find nothing at any level and the two pieces are genuinely disconnected.
//
// LOWER BOUND. reading the ops is Omega(q); you can't answer online without seeing
// each op. HDLT spends O(log^2 n) amortized per update and O(log n) per query (the
// treap depth of F_0). the offline segment-tree-on-time method beats it at
// O(log q * alpha) per op -- but only because it gets the future for free. online,
// the level machinery is the price of that missing future, and it's known to be
// near-optimal: any structure has an Omega(log n / log log n) amortized cell-probe
// bound per op (Patrascu-Demaine). so we're a log factor off the theoretical floor
// and there's no segment tree to rescue us.

// ---------------------------------------------------------------------------
// EulerTourForest -- a forest of Euler-tour treaps over a fixed vertex set. link,
// cut, connectivity, component size, and two families of aggregate bits used by
// the level search. this is the layer to trust before the levels go on top: it has
// a clean link/cut/connected contract you can differential-test on its own.
// ---------------------------------------------------------------------------
class EulerTourForest {
public:
    // an edge handle -- the two arc tokens a link created. cut takes it back.
    struct Edge {
        int a = 0;  // token for direction u -> v
        int b = 0;  // token for direction v -> u
    };

    // n vertices, 1-based. every vertex starts isolated: its tour is the single
    // vertex node [V_u]. node 0 is the null sentinel -- size 0, all aggregates off.
    void init(int n) {
        n_ = n;
        std::size_t cap = static_cast<std::size_t>(n) + 1;
        l_.assign(cap, 0);
        r_.assign(cap, 0);
        p_.assign(cap, 0);
        sz_.assign(cap, 0);
        vcnt_.assign(cap, 0);
        pr_.assign(cap, 0);
        is_vtx_.assign(cap, 0);
        vtx_.assign(cap, 0);
        edge_.assign(cap, 0);
        mv_.assign(cap, 0);
        avm_.assign(cap, 0);
        me_.assign(cap, 0);
        aem_.assign(cap, 0);
        free_.clear();
        rng_ = 0x9E3779B97F4A7C15ULL;
        for (int u = 1; u <= n; ++u) {
            std::size_t x = static_cast<std::size_t>(u);
            is_vtx_[x] = 1;
            vtx_[x] = u;
            sz_[x] = 1;
            vcnt_[x] = 1;
            pr_[x] = next_pri();
        }
    }

    // the treap root of u's component -- the component identity token.
    int root_of(int u) const { return find_root(u); }

    bool connected(int u, int v) const {
        if (u == v) return true;
        return find_root(u) == find_root(v);
    }

    // vertices in u's component. a tour of a k-vertex tree holds exactly k vertex
    // nodes, so the aggregate count reads size straight off the root.
    int comp_size(int u) const {
        return vcnt_[static_cast<std::size_t>(find_root(u))];
    }

    // link u and v -- they live in different components by contract. reroot each
    // tour to open at its endpoint, then stitch: tour_u, arc(u->v), tour_v,
    // arc(v->u) is a closed Euler walk of the joined tree. edge_id is the caller's
    // handle for the edge, stamped on both tokens so find_marked_edge hands it back.
    Edge link(int u, int v, int edge_id) {
        int su = reroot_at(u);
        int sv = reroot_at(v);
        int a = alloc_token(edge_id);
        int b = alloc_token(edge_id);
        merge(merge(merge(su, a), sv), b);  // root is found on demand.
        return {a, b};
    }

    // cut the edge whose tokens are e.a, e.b. the tour is cyclic, so first ROTATE
    // it to open at a -- merge(after-a, before-a) -- which wraps b forward past a no
    // matter which came first linearly. then peel a, split at b: the block between
    // them is v's side of the tree, everything after b is u's. both tokens freed.
    //
    // the rotation is the crux. without it, when b precedes a in the linear order,
    // b stays stranded in the before-a half and the forward split at b is nonsense --
    // it leaves a freed token wired in as a live child and rots the aggregates.
    void cut(Edge e) {
        int t = find_root(e.a);
        int ka = order(e.a);
        auto [p1, p2] = split(t, ka);         // p2 opens at a, p1 = everything before.
        int rot = merge(p2, p1);              // cyclic rotate: tour now opens at a.
        auto [aonly, rest] = split(rot, 1);   // aonly == {a}.
        int kb = order(e.b);                  // b's position inside rest (forward from a).
        auto [mid, rest2] = split(rest, kb);  // mid == v-side sub-tour.
        auto [bonly, z] = split(rest2, 1);    // bonly == {b}; z == u-side tour.
        (void)aonly;
        (void)mid;
        (void)bonly;
        (void)z;
        free_token(e.a);
        free_token(e.b);
    }

    // ---- aggregate bits -----------------------------------------------------
    // mv: this vertex has an incident non-tree edge at the forest's level. set it
    // on the stable vertex node, then refresh the aggregate up to the root.
    void set_vertex_mark(int u, bool on) {
        std::size_t x = static_cast<std::size_t>(u);
        char c = on ? 1 : 0;
        if (mv_[x] == c) return;
        mv_[x] = c;
        refresh_up(u);
    }
    bool vertex_mark(int u) const {
        return mv_[static_cast<std::size_t>(u)] != 0;
    }

    // me: this arc token's edge is a tree edge OF the forest's level (push-up
    // candidate). set on a token, refreshed up.
    void set_edge_mark(int token, bool on) {
        std::size_t x = static_cast<std::size_t>(token);
        char c = on ? 1 : 0;
        if (me_[x] == c) return;
        me_[x] = c;
        refresh_up(token);
    }

    // any vertex in u's component carrying mv, or 0. descends by the avm aggregate
    // in O(log n).
    int find_marked_vertex(int u) const {
        int x = find_root(u);
        if (x == 0 || !avm_[static_cast<std::size_t>(x)]) return 0;
        for (;;) {
            std::size_t xs = static_cast<std::size_t>(x);
            int lc = l_[xs];
            if (lc && avm_[static_cast<std::size_t>(lc)]) {
                x = lc;
                continue;
            }
            if (is_vtx_[xs] && mv_[xs]) return vtx_[xs];
            x = r_[xs];  // aggregate guarantees the mark is to the right.
        }
    }

    // any tree edge id in u's component whose token carries me, or -1. descends by
    // the aem aggregate. returns the edge id stored on the token.
    int find_marked_edge(int u) const {
        int x = find_root(u);
        if (x == 0 || !aem_[static_cast<std::size_t>(x)]) return -1;
        for (;;) {
            std::size_t xs = static_cast<std::size_t>(x);
            int lc = l_[xs];
            if (lc && aem_[static_cast<std::size_t>(lc)]) {
                x = lc;
                continue;
            }
            if (me_[xs]) return edge_[xs];
            x = r_[xs];
        }
    }

private:
    int n_ = 0;

    std::vector<int> l_, r_, p_;
    std::vector<int> sz_, vcnt_;
    std::vector<std::uint32_t> pr_;
    std::vector<char> is_vtx_;
    std::vector<int> vtx_;
    std::vector<int> edge_;  // token -> the edge id it belongs to.
    std::vector<char> mv_, avm_, me_, aem_;
    std::vector<int> free_;  // freed token slots, reused before growing.
    std::uint64_t rng_ = 0;

    std::uint32_t next_pri() {
        // splitmix64 -- fast, deterministic, good enough spread for treap balance.
        rng_ += 0x9E3779B97F4A7C15ULL;
        std::uint64_t z = rng_;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        return static_cast<std::uint32_t>(z);
    }

    int alloc_token(int edge_id) {
        int x;
        if (!free_.empty()) {
            x = free_.back();
            free_.pop_back();
        } else {
            x = static_cast<int>(l_.size());
            l_.push_back(0);
            r_.push_back(0);
            p_.push_back(0);
            sz_.push_back(0);
            vcnt_.push_back(0);
            pr_.push_back(0);
            is_vtx_.push_back(0);
            vtx_.push_back(0);
            edge_.push_back(0);
            mv_.push_back(0);
            avm_.push_back(0);
            me_.push_back(0);
            aem_.push_back(0);
        }
        std::size_t xs = static_cast<std::size_t>(x);
        l_[xs] = r_[xs] = p_[xs] = 0;
        sz_[xs] = 1;
        vcnt_[xs] = 0;  // a token is not a vertex -- it doesn't count toward size.
        pr_[xs] = next_pri();
        is_vtx_[xs] = 0;
        vtx_[xs] = 0;
        edge_[xs] = edge_id;
        mv_[xs] = avm_[xs] = me_[xs] = aem_[xs] = 0;
        return x;
    }

    void free_token(int x) {
        std::size_t xs = static_cast<std::size_t>(x);
        l_[xs] = r_[xs] = p_[xs] = 0;
        me_[xs] = aem_[xs] = 0;
        free_.push_back(x);
    }

    int size_of(int x) const { return sz_[static_cast<std::size_t>(x)]; }

    // recombine a node's aggregates from its two children. size, vertex count, and
    // the two mark ORs -- everything the descents read.
    void pull(int x) {
        if (x == 0) return;
        std::size_t xs = static_cast<std::size_t>(x);
        int lc = l_[xs], rc = r_[xs];
        std::size_t ls = static_cast<std::size_t>(lc);
        std::size_t rs = static_cast<std::size_t>(rc);
        sz_[xs] = 1 + sz_[ls] + sz_[rs];
        vcnt_[xs] = (is_vtx_[xs] ? 1 : 0) + vcnt_[ls] + vcnt_[rs];
        avm_[xs] =
            static_cast<char>((is_vtx_[xs] && mv_[xs]) || avm_[ls] || avm_[rs]);
        aem_[xs] = static_cast<char>(me_[xs] || aem_[ls] || aem_[rs]);
    }

    int find_root(int x) const {
        while (p_[static_cast<std::size_t>(x)]) x = p_[static_cast<std::size_t>(x)];
        return x;
    }

    // 0-indexed position of x within its treap, by summing left-subtree sizes on
    // the walk to the root.
    int order(int x) const {
        int k = size_of(l_[static_cast<std::size_t>(x)]);
        while (p_[static_cast<std::size_t>(x)]) {
            int par = p_[static_cast<std::size_t>(x)];
            if (r_[static_cast<std::size_t>(par)] == x)
                k += size_of(l_[static_cast<std::size_t>(par)]) + 1;
            x = par;
        }
        return k;
    }

    // split t into the first k nodes and the rest. both returned roots have a null
    // parent. parents of moved children are kept consistent.
    std::pair<int, int> split(int t, int k) {
        if (t == 0) return {0, 0};
        std::size_t ts = static_cast<std::size_t>(t);
        int ls = size_of(l_[ts]);
        if (ls >= k) {
            auto [a, b] = split(l_[ts], k);
            l_[ts] = b;
            if (b) p_[static_cast<std::size_t>(b)] = t;
            pull(t);
            p_[ts] = 0;
            if (a) p_[static_cast<std::size_t>(a)] = 0;
            return {a, t};
        } else {
            auto [a, b] = split(r_[ts], k - ls - 1);
            r_[ts] = a;
            if (a) p_[static_cast<std::size_t>(a)] = t;
            pull(t);
            p_[ts] = 0;
            if (b) p_[static_cast<std::size_t>(b)] = 0;
            return {t, b};
        }
    }

    int merge(int a, int b) {
        if (!a) return b;
        if (!b) return a;
        if (pr_[static_cast<std::size_t>(a)] < pr_[static_cast<std::size_t>(b)]) {
            int m = merge(r_[static_cast<std::size_t>(a)], b);
            r_[static_cast<std::size_t>(a)] = m;
            if (m) p_[static_cast<std::size_t>(m)] = a;
            pull(a);
            p_[static_cast<std::size_t>(a)] = 0;
            return a;
        } else {
            int m = merge(a, l_[static_cast<std::size_t>(b)]);
            l_[static_cast<std::size_t>(b)] = m;
            if (m) p_[static_cast<std::size_t>(m)] = b;
            pull(b);
            p_[static_cast<std::size_t>(b)] = 0;
            return b;
        }
    }

    // rotate u's tour to open at u's vertex node, and return the new root.
    int reroot_at(int u) {
        int t = find_root(u);
        int k = order(u);
        auto [pre, suf] = split(t, k);  // suf opens at u.
        return merge(suf, pre);
    }

    // a mark on x changed -- recompute aggregates from x up to the root.
    void refresh_up(int x) {
        while (x) {
            pull(x);
            x = p_[static_cast<std::size_t>(x)];
        }
    }
};

// ---------------------------------------------------------------------------
// DynamicConnectivity -- HDLT on top of L+1 Euler-tour forests.
// ---------------------------------------------------------------------------
class DynamicConnectivity {
public:
    // n vertices, 1-based. sizes the level structure and the L+1 forests.
    void init(int n) {
        n_ = n;
        levels_ = 1;
        while ((1 << levels_) < n) ++levels_;  // L = ceil(log2 n), forests 0..L.
        forest_.assign(static_cast<std::size_t>(levels_ + 1), EulerTourForest());
        for (auto& f : forest_) f.init(n);

        adj_.assign(static_cast<std::size_t>(levels_ + 1), {});
        for (auto& lvl : adj_)
            lvl.assign(static_cast<std::size_t>(n) + 1, {});

        edges_.clear();
        edges_.push_back(EdgeRec{});  // id 0 unused.
        edge_id_.clear();
    }

    void insert(int u, int v) {
        int id = static_cast<int>(edges_.size());
        EdgeRec rec;
        rec.u = u;
        rec.v = v;
        edges_.push_back(rec);
        edge_id_[key(u, v)] = id;

        if (!forest_[0].connected(u, v)) {
            make_tree(id, 0);   // spans two components -- a level-0 tree edge.
        } else {
            add_nontree(id, 0);  // already connected -- a level-0 non-tree edge.
        }
    }

    void remove(int u, int v) {
        auto it = edge_id_.find(key(u, v));
        int id = it->second;
        edge_id_.erase(it);
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];

        if (!e.tree) {
            drop_nontree(id);
            return;
        }

        int lv = e.level;
        // pull the tree edge out of every forest that holds it: F_0..F_lv.
        for (int i = 0; i <= lv; ++i)
            forest_[static_cast<std::size_t>(i)].cut(e.tok[i]);
        e.tree = false;

        // hunt a replacement from the edge's level down to 0.
        for (int i = lv; i >= 0; --i)
            if (try_reconnect(i, u, v)) break;
    }

    bool connected(int u, int v) const { return forest_[0].connected(u, v); }

private:
    struct EdgeRec {
        int u = 0;
        int v = 0;
        int level = 0;
        bool tree = false;
        EulerTourForest::Edge tok[32];  // per-level tokens while a tree edge.
        int pos_u = -1;                 // index in adj_[level][u] as a non-tree edge.
        int pos_v = -1;                 // index in adj_[level][v].
    };

    int n_ = 0;
    int levels_ = 0;  // L. forests indexed 0..L.
    std::vector<EulerTourForest> forest_;
    std::vector<std::vector<std::vector<int>>> adj_;  // adj_[level][vertex] -> ids.
    std::vector<EdgeRec> edges_;

    struct KeyHash {
        std::size_t operator()(std::uint64_t k) const {
            k ^= k >> 33;
            k *= 0xFF51AFD7ED558CCDULL;
            k ^= k >> 33;
            return static_cast<std::size_t>(k);
        }
    };
    std::unordered_map<std::uint64_t, int, KeyHash> edge_id_;  // present edge -> id.

    static std::uint64_t key(int u, int v) {
        int a = u < v ? u : v;
        int b = u < v ? v : u;
        return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
    }

    // make edge id a tree edge at level lv: splice it into F_0..F_lv and stamp the
    // "tree edge of level lv" mark on its F_lv tokens so the push-up scan finds it.
    void make_tree(int id, int lv) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        e.tree = true;
        e.level = lv;
        for (int i = 0; i <= lv; ++i)
            e.tok[i] = forest_[static_cast<std::size_t>(i)].link(e.u, e.v, id);
        mark_tree_level(id, lv, true);
    }

    // record a non-tree edge at level lv in both incidence lists; light the mv bit
    // on each endpoint in F_lv.
    void add_nontree(int id, int lv) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        e.tree = false;
        e.level = lv;
        auto& au = adj_[static_cast<std::size_t>(lv)][static_cast<std::size_t>(e.u)];
        auto& av = adj_[static_cast<std::size_t>(lv)][static_cast<std::size_t>(e.v)];
        e.pos_u = static_cast<int>(au.size());
        au.push_back(id);
        e.pos_v = static_cast<int>(av.size());
        av.push_back(id);
        forest_[static_cast<std::size_t>(lv)].set_vertex_mark(e.u, true);
        forest_[static_cast<std::size_t>(lv)].set_vertex_mark(e.v, true);
    }

    // remove a non-tree edge from its incidence lists (swap-pop), clearing the mv
    // bit on any endpoint whose list just emptied.
    void drop_nontree(int id) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        int lv = e.level;
        remove_from_adj(lv, e.u, e.pos_u, id);
        remove_from_adj(lv, e.v, e.pos_v, id);
        maybe_clear_mark(lv, e.u);
        maybe_clear_mark(lv, e.v);
    }

    void remove_from_adj(int lv, int vtx, int pos, int id) {
        auto& lst = adj_[static_cast<std::size_t>(lv)][static_cast<std::size_t>(vtx)];
        int last = static_cast<int>(lst.size()) - 1;
        int moved = lst[static_cast<std::size_t>(last)];
        lst[static_cast<std::size_t>(pos)] = moved;
        lst.pop_back();
        if (moved != id) {
            // fix the moved edge's stored position for this endpoint.
            EdgeRec& m = edges_[static_cast<std::size_t>(moved)];
            if (m.u == vtx && m.pos_u == last)
                m.pos_u = pos;
            else
                m.pos_v = pos;
        }
    }

    void maybe_clear_mark(int lv, int vtx) {
        if (adj_[static_cast<std::size_t>(lv)][static_cast<std::size_t>(vtx)].empty())
            forest_[static_cast<std::size_t>(lv)].set_vertex_mark(vtx, false);
    }

    // stamp / clear the "tree edge of level lv" mark on an edge's F_lv tokens.
    void mark_tree_level(int id, int lv, bool on) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        forest_[static_cast<std::size_t>(lv)].set_edge_mark(e.tok[lv].a, on);
        forest_[static_cast<std::size_t>(lv)].set_edge_mark(e.tok[lv].b, on);
    }

    // try to reconnect the two pieces at level i. returns true if a replacement was
    // installed, false to fall through to level i-1.
    bool try_reconnect(int i, int u, int v) {
        EulerTourForest& F = forest_[static_cast<std::size_t>(i)];
        int su = F.comp_size(u);
        int sv = F.comp_size(v);
        int small_v = (su <= sv) ? u : v;       // a vertex in the small piece.
        int other_root = (su <= sv) ? F.root_of(v) : F.root_of(u);  // far root.

        // 1. push every level-i tree edge of the small piece up to level i+1.
        for (;;) {
            int eid = F.find_marked_edge(small_v);
            if (eid < 0) break;
            promote_tree_edge(eid, i);
        }

        // 2. scan the small piece's level-i non-tree edges for a crossing one.
        for (;;) {
            int w = F.find_marked_vertex(small_v);
            if (w == 0) break;
            auto& lst = adj_[static_cast<std::size_t>(i)][static_cast<std::size_t>(w)];
            while (!lst.empty()) {
                int eid = lst.back();
                EdgeRec& e = edges_[static_cast<std::size_t>(eid)];
                int x = (e.u == w) ? e.v : e.u;  // the far endpoint.
                if (F.root_of(x) == other_root) {
                    // crosses the cut -- the replacement. promote to a tree edge at
                    // level i, spanning F_0..F_i.
                    drop_nontree(eid);
                    make_tree(eid, i);
                    return true;
                }
                raise_nontree(eid, i);  // stays inside -- push it up, keep scanning.
            }
        }
        return false;
    }

    // raise a level-i tree edge to level i+1: add it to F_{i+1}, move its level mark
    // from F_i up to F_{i+1}.
    void promote_tree_edge(int id, int i) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        mark_tree_level(id, i, false);
        e.tok[i + 1] = forest_[static_cast<std::size_t>(i + 1)].link(e.u, e.v, id);
        e.level = i + 1;
        mark_tree_level(id, i + 1, true);
    }

    // raise a level-i non-tree edge to level i+1: move it between incidence lists
    // and refresh the mv bits in both forests.
    void raise_nontree(int id, int i) {
        EdgeRec& e = edges_[static_cast<std::size_t>(id)];
        int u = e.u, v = e.v;
        remove_from_adj(i, u, e.pos_u, id);
        remove_from_adj(i, v, e.pos_v, id);
        maybe_clear_mark(i, u);
        maybe_clear_mark(i, v);
        e.level = i + 1;
        auto& au = adj_[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(u)];
        auto& av = adj_[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(v)];
        e.pos_u = static_cast<int>(au.size());
        au.push_back(id);
        e.pos_v = static_cast<int>(av.size());
        av.push_back(id);
        forest_[static_cast<std::size_t>(i + 1)].set_vertex_mark(u, true);
        forest_[static_cast<std::size_t>(i + 1)].set_vertex_mark(v, true);
    }
};

}  // namespace p0046
