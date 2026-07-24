#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace p0058 {

// the top tree -- the hardest dynamic-tree structure there is, a full tier above
// the link-cut tree. a forest, fully dynamic: link two trees, cut an edge, and ask
// for the DIAMETER of a whole tree. all three in O(log n) amortized.
//
// WHY A TOP TREE, NOT AN LCT. an LCT exposes a root-to-node PATH and folds it. the
// diameter is not a path aggregate -- it's the longest path between ANY two vertices,
// and its endpoints can sit in two different side-branches that never share a
// root-to-node path. so a structure that can only see one path at a time can't read
// it. you need to fold the ENTIRE tree, branches and all, into one root cluster.
// that is exactly what a top tree does, and rake -- not just compress -- is the
// operation that folds the branches in.
//
// THE CLUSTER. a cluster is a connected piece of the tree whose interface to the
// rest is at most two BOUNDARY vertices. every cluster stores, over the piece it
// covers:
//   len  -- the weighted length of the path between its two boundary vertices,
//   d1   -- max distance from boundary 1 to any vertex inside,
//   d2   -- max distance from boundary 2 to any vertex inside,
//   diam -- the longest path lying fully inside.
// a single edge {u, v, w} is the base cluster: boundaries u, v; len = d1 = d2 =
// diam = w (from either boundary the far vertex is w away, and the only path is w).
//
// COMPRESS -- glue two path clusters x (a..m) and y (m..b) at a shared degree-2
// boundary m into one path cluster a..b:
//   len  = len_x + len_y
//   d1   = max(d1_x, len_x + d1_y)              // from a: stay in x, or cross to y
//   d2   = max(d2_y, len_y + d2_x)              // from b: stay in y, or cross to x
//   diam = max(diam_x, diam_y, d2_x + d1_y)     // or a path that bends through m
// the d2_x + d1_y term is the whole point: the longest path through m is the deepest
// reach into x from m plus the deepest reach into y from m.
//
// RAKE -- a subtree hangs off a boundary. fold that hanging cluster c into a path
// cluster p, keeping p's two boundaries. c is oriented so its interface vertex is
// on the "right" (its d2 faces the shared boundary). raking c onto p's right
// boundary:
//   len  = len_p                               // the path between p's boundaries is unchanged
//   d1   = max(d1_p, len_p + d2_c)             // from p's left: cross p, then dive into c
//   d2   = max(d2_p, d2_c)                     // from the shared boundary: into p, or into c
//   diam = max(diam_p, diam_c, d2_p + d2_c)    // or a path that bends through the boundary
// the d2_p + d2_c term is the same idea as compress's bend, but the two deep reaches
// now both hang off the one shared boundary instead of a through-path.
//
// so a query is: fold the tree to its root cluster, read diam. that fold is the
// forest-oracle's O(component) BFS+DFS collapsed into the tree's O(log n) height.
//
// THE STRUCTURE -- a self-adjusting top tree (Tarjan-Werneck style). two kinds of
// splay trees interleave. COMPRESS trees: a splay tree whose in-order is a maximal
// path of the contracted tree; every internal node compresses its two children.
// RAKE trees: a splay tree hung off a compress node, collecting every subtree raked
// onto that node's middle vertex. base clusters -- edges -- are the leaves of both.
// expose(v) rotates and re-hangs O(log n) clusters so v lands on the root cluster's
// boundary; link and cut splice a constant number of clusters and re-fold. the
// amortized O(log n) is the top-tree access lemma, the same shape that pays for the
// LCT's preferred-path splays.
//
// COMPLEXITY. O((n + q) log n) total. reading the ops alone is Omega(q), so the log
// factor is the only slack. the naive alternative -- BFS the component and DFS its
// diameter per query -- is O(q * n); on 1e5 vertices with 2e5 queries that is 2e10
// versus the top tree's ~2e5 * 17. the win is the height of the tree.
//
// INT64. one path holds at most n-1 edges, each weight <= 1e9, so a diameter is at
// most ~1e5 * 1e9 = 1e14 -- inside int64's 9.2e18 with four orders of magnitude to
// spare. every len/d1/d2/diam is a sum of edge weights along a real path, so the
// same bound covers all of them.
//
// INVARIANTS (what fix() restores after every splice, and what the whole thing
// leans on):
//   I1. a compress node's two boundaries are its left child's boundary 0 and its
//       right child's boundary 1; the children meet at one shared vertex
//       (left.boundary1 == right.boundary0) -- the compressed degree-2 vertex.
//   I2. a compress node's rake tree holds exactly the subtrees hanging off that
//       shared middle vertex; fold(node) rakes the rake-tree into the left child,
//       then compresses with the right child.
//   I3. a rake node's boundaries are its left child's boundaries; its fold rakes
//       the two children (both hanging off the same vertex).
//   I4. handle[v] points to the one cluster that "owns" v -- the compress node whose
//       middle vertex is v, or, for a boundary vertex, the root cluster of its tree.
//       expose starts from handle[v] and climbs.
//   I5. a lazy reverse flag on a compress node swaps its two boundaries and its
//       d1/d2; it commutes with fold because diam and len are reversal-symmetric and
//       d1/d2 simply swap. push() ships it to the children before any structural read.
template <class Cl>
class TopForest {
public:
    using F = typename Cl::F;

    // n real vertices, 1..n externally. each real vertex is pinned to a private
    // dummy vertex by a zero-length edge -- a leaf that hangs off it and never moves.
    // that keeps every real vertex's handle non-null at all times, so link and cut
    // never special-case an isolated point, and a zero edge can't lengthen any path
    // or inflate a diameter. cost is one extra vertex and edge per real vertex.
    explicit TopForest(int n) : n_(n) {
        std::size_t verts = static_cast<std::size_t>(2 * n);
        handle_.assign(verts, -1);
        // rough upper bound on nodes ever created: each of the n dummy links and each
        // future op adds O(1) clusters. reserve generously to dodge reallocation.
        nd_.reserve(verts * 4 + 64);
        for (int i = 0; i < n; ++i)
            makeEdge(i, n + i, Cl::id());  // real i -- dummy (n+i), weight 0.
    }

    // add edge {u, v} of weight w. u, v are 1-based and in different trees.
    void link(int u, int v, std::int64_t w) {
        linkRaw(u - 1, v - 1, Cl::base(w));
    }

    // remove edge {u, v}. it exists by contract. 1-based.
    void cut(int u, int v) { cutRaw(u - 1, v - 1); }

    // whole-tree fold for the tree containing 1-based u. read .diam off it for the
    // diameter. expose(u) collapses u's whole tree into one root cluster.
    F rootAgg(int u) { return nd_[static_cast<std::size_t>(expose(u - 1))].fold; }

    // fold of the path between 1-based u and v (they must be connected). used to
    // gate the path-length layer against its own oracle before diameter goes on top.
    F pathAgg(int u, int v) { return pathAggRaw(u - 1, v - 1); }

private:
    enum Kind : std::uint8_t { EDGE = 0, COMPRESS = 1, RAKE = 2 };

    // one node covers all three roles -- edge leaf, compress, rake -- with uniform
    // fields so endpoints() and fold read the same way everywhere. ch are compress
    // children (COMPRESS) or rake children (RAKE); rake is the hung rake-tree root
    // (COMPRESS only). v0/v1 are the two boundary vertices. par is the parent cluster.
    struct Nd {
        Kind kind;
        bool rev = false;    // lazy path reverse (COMPRESS).
        bool guard = false;  // soft-expose fence: hide this node from parent_dir.
        int ch[2] = {-1, -1};
        int rake = -1;
        int par = -1;
        int v0 = -1, v1 = -1;
        F fold;
    };

    int n_;
    std::vector<Nd> nd_;      // the cluster arena. index == node id; -1 is null.
    std::vector<int> handle_; // handle_[vertex] -> owning cluster (invariant I4).

    // ---- arena + uniform accessors -----------------------------------------
    int newNode(Kind k) {
        int id = static_cast<int>(nd_.size());
        nd_.emplace_back();
        nd_[static_cast<std::size_t>(id)].kind = k;
        return id;
    }
    Nd& N(int id) { return nd_[static_cast<std::size_t>(id)]; }
    // boundary vertex on side dir (0 or 1) of any node.
    int endp(int id, int dir) { return dir == 0 ? N(id).v0 : N(id).v1; }
    void setHandle(int v, int id) { handle_[static_cast<std::size_t>(v)] = id; }

    // a fresh edge leaf {a, b, weight}. its own fix pins the handles.
    int makeEdge(int a, int b, F weight) {
        int e = newNode(EDGE);
        N(e).v0 = a;
        N(e).v1 = b;
        N(e).fold = weight;
        fix(e);
        return e;
    }
    int makeCompress(int l, int r) {
        int c = newNode(COMPRESS);
        N(c).ch[0] = l;
        N(c).ch[1] = r;
        N(c).fold = Cl::id();
        fix(c);
        return c;
    }
    int makeRake(int l, int r) {
        int c = newNode(RAKE);
        N(c).ch[0] = l;
        N(c).ch[1] = r;
        N(c).fold = Cl::id();
        fix(c);
        return c;
    }

    // ---- fix / push / reverse, dispatched by kind --------------------------
    void fix(int id) {
        switch (N(id).kind) {
            case EDGE: fixEdge(id); break;
            case COMPRESS: fixComp(id); break;
            case RAKE: fixRake(id); break;
        }
    }
    void pushNode(int id) {
        if (N(id).kind == COMPRESS) pushComp(id);
        // edge and rake have nothing pending.
    }
    void reverseNode(int id) {
        Nd& t = N(id);
        switch (t.kind) {
            case EDGE:
                std::swap(t.v0, t.v1);
                Cl::rev(t.fold);
                break;
            case COMPRESS:
                std::swap(t.v0, t.v1);
                Cl::rev(t.fold);
                t.rev = !t.rev;
                break;
            case RAKE: break;  // a rake node has no orientation to flip.
        }
    }

    void fixEdge(int id) {
        Nd& t = N(id);
        int p = t.par;
        if (p != -1 && N(p).kind == COMPRESS) {
            if (!parentDirComp(id).ok) setHandle(t.v0, id);
        } else if (p != -1 && N(p).kind == RAKE) {
            setHandle(t.v0, id);
        } else {  // a lone root edge owns both its boundaries.
            setHandle(t.v0, id);
            setHandle(t.v1, id);
        }
    }

    void fixComp(int id) {
        pushComp(id);
        Nd& t = N(id);
        int l = t.ch[0], r = t.ch[1];
        t.v0 = endp(l, 0);
        t.v1 = endp(r, 1);
        // fold: rake the hung subtrees into the left child at the middle vertex,
        // then compress with the right child (invariant I2).
        F left = (t.rake != -1) ? Cl::rake(N(l).fold, N(t.rake).fold) : N(l).fold;
        t.fold = Cl::comp(left, N(r).fold);
        setHandle(endp(l, 1), id);  // the compressed middle vertex (I1).
        int p = t.par;
        if (p != -1 && N(p).kind == COMPRESS) {
            if (!parentDirComp(id).ok) setHandle(t.v0, id);
        } else if (p != -1 && N(p).kind == RAKE) {
            setHandle(t.v0, id);
        } else {
            setHandle(t.v0, id);
            setHandle(t.v1, id);
        }
    }

    void fixRake(int id) {
        Nd& t = N(id);
        int l = t.ch[0], r = t.ch[1];
        t.v0 = endp(l, 0);
        t.v1 = endp(l, 1);  // a rake node inherits the retained child's boundaries (I3).
        t.fold = Cl::rake(N(l).fold, N(r).fold);
    }

    void pushComp(int id) {
        Nd& t = N(id);
        if (t.rev) {
            std::swap(t.ch[0], t.ch[1]);
            reverseNode(t.ch[0]);
            reverseNode(t.ch[1]);
            t.rev = false;
        }
    }

    // ---- parent-direction probes -------------------------------------------
    // the boundary between splay trees is drawn by these: a node is a genuine child
    // of its parent only if the parent lists it in ch[]. a guard node lies about its
    // children (returns none) so soft_expose can fence a subtree off mid-operation.
    struct PD {
        bool ok;
        int dir;
        int p;
    };
    PD parentDirComp(int child) {  // child is a compress/edge, parent a non-guard compress.
        int p = N(child).par;
        if (p != -1 && N(p).kind == COMPRESS && !N(p).guard) {
            if (N(p).ch[0] == child) return {true, 0, p};
            if (N(p).ch[1] == child) return {true, 1, p};
        }
        return {false, 0, -1};
    }
    PD parentDirCompGuard(int child) {  // same, but a guard still owns its children.
        int p = N(child).par;
        if (p != -1 && N(p).kind == COMPRESS) {
            if (N(p).ch[0] == child) return {true, 0, p};
            if (N(p).ch[1] == child) return {true, 1, p};
        }
        return {false, 0, -1};
    }
    PD parentDirCompRake(int child) {  // child is a comp node sitting as a rake leaf.
        int p = N(child).par;
        if (p != -1 && N(p).kind == RAKE) {
            if (N(p).ch[0] == child) return {true, 0, p};
            if (N(p).ch[1] == child) return {true, 1, p};
        }
        return {false, 0, -1};
    }
    PD parentDirRake(int child) {  // child is a rake node, parent a rake node.
        int p = N(child).par;
        if (p != -1 && N(p).kind == RAKE) {
            if (N(p).ch[0] == child) return {true, 0, p};
            if (N(p).ch[1] == child) return {true, 1, p};
        }
        return {false, 0, -1};
    }

    // ---- compress-tree splay -----------------------------------------------
    // t is child(dir) of x; rotate t above x. the four tails re-hook t into whatever
    // held x -- a compress parent, a rake parent (t becomes a rake leaf), the rake
    // pointer of a compress (t becomes its rake root), or nothing (t is a new root).
    void rotateComp(int t, int x, int dir) {
        int y = N(x).par;
        PD par = parentDirCompGuard(x);
        PD rpar = parentDirCompRake(x);

        int c = N(t).ch[dir];
        N(x).ch[dir ^ 1] = c;
        N(c).par = x;
        N(t).ch[dir] = x;
        N(x).par = t;
        N(t).par = y;

        if (par.ok) {
            N(par.p).ch[par.dir] = t;
            fix(x);
            fix(t);
            if (!N(par.p).guard) fix(par.p);
        } else if (rpar.ok) {
            N(rpar.p).ch[rpar.dir] = t;
            fix(x);
            fix(t);
            fix(rpar.p);
        } else if (y != -1 && N(y).kind == COMPRESS) {
            N(y).rake = t;  // x was y's rake root; t inherits that slot.
            fix(x);
            fix(t);
            if (!N(y).guard) fix(y);
        } else {
            fix(x);
            fix(t);
        }
    }

    void splayComp(int t) {
        pushComp(t);
        for (;;) {
            PD q = parentDirComp(t);
            if (!q.ok) break;
            PD r = parentDirComp(q.p);
            if (r.ok) {
                int rp = N(r.p).par;
                if (rp != -1) pushNode(rp);
                pushNode(r.p);
                pushNode(q.p);
                pushComp(t);
                // recompute directions AFTER the pushes -- a reverse can swap them.
                int qt = parentDirComp(t).dir;
                int rq = parentDirComp(q.p).dir;
                if (rq == qt) {
                    rotateComp(q.p, r.p, rq ^ 1);
                    rotateComp(t, q.p, qt ^ 1);
                } else {
                    rotateComp(t, q.p, qt ^ 1);
                    rotateComp(t, r.p, rq ^ 1);
                }
            } else {
                int qp = N(q.p).par;
                if (qp != -1) pushNode(qp);
                pushComp(q.p);
                pushComp(t);
                int qt = parentDirComp(t).dir;
                rotateComp(t, q.p, qt ^ 1);
            }
        }
    }

    // ---- rake-tree splay ---------------------------------------------------
    void rotateRake(int t, int x, int dir) {
        int y = N(x).par;
        PD par = parentDirRake(x);

        int c = N(t).ch[dir];
        N(x).ch[dir ^ 1] = c;
        N(c).par = x;
        N(t).ch[dir] = x;
        N(x).par = t;
        N(t).par = y;

        if (par.ok) {
            N(par.p).ch[par.dir] = t;
            fix(x);
            fix(t);
            fix(par.p);
        } else if (y != -1 && N(y).kind == COMPRESS) {
            N(y).rake = t;  // x was a compress's rake root.
            fix(x);
            fix(t);
            if (!N(y).guard) fix(y);
        } else {
            fix(x);
            fix(t);
        }
    }

    void splayRake(int t) {
        for (;;) {
            PD q = parentDirRake(t);
            if (!q.ok) break;
            PD r = parentDirRake(q.p);
            if (r.ok) {
                int qt = parentDirRake(t).dir;
                int rq = parentDirRake(q.p).dir;
                if (rq == qt) {
                    rotateRake(q.p, r.p, rq ^ 1);
                    rotateRake(t, q.p, qt ^ 1);
                } else {
                    rotateRake(t, q.p, qt ^ 1);
                    rotateRake(t, r.p, rq ^ 1);
                }
            } else {
                int qt = parentDirRake(t).dir;
                rotateRake(t, q.p, qt ^ 1);
            }
        }
    }

    // ---- expose ------------------------------------------------------------
    // lift `node` to the root cluster of its tree, climbing splay tree to splay tree.
    // each turn: splay within the current compress tree, walk up through a rake tree
    // (splay it, its parent is the compress above) or a compress parent, then swap the
    // path we came from into that compress as a child, demoting the child it displaces
    // into the rake tree. that demote-into-rake is where branches get folded in.
    int exposeRaw(int node) {
        for (;;) {
            if (N(node).kind == COMPRESS) splayComp(node);
            int p = N(node).par;
            if (p == -1) break;

            int n;
            if (N(p).kind == RAKE) {
                pushNode(p);
                splayRake(p);
                n = N(p).par;  // a rake tree's root hangs off a compress (I2).
            } else {  // COMPRESS
                pushComp(p);
                if (N(p).guard && parentDirCompGuard(node).ok) break;  // fenced by soft_expose.
                n = p;
            }

            splayComp(n);
            PD nd = parentDirCompGuard(n);
            int dir = nd.ok ? nd.dir : 0;
            if (dir == 1) {
                reverseNode(N(n).ch[dir]);
                pushNode(N(n).ch[dir]);
                reverseNode(node);
                pushNode(node);
            }

            PD rp = parentDirCompRake(node);  // was `node` already a rake leaf under n?
            if (rp.ok) {
                int nch = N(n).ch[dir];
                pushNode(nch);
                pushNode(rp.p);
                N(rp.p).ch[rp.dir] = nch;  // demoted child takes node's old rake slot.
                N(nch).par = rp.p;
                N(n).ch[dir] = node;  // node is promoted into the compress path.
                N(node).par = n;
                fix(nch);
                fix(rp.p);
                fix(node);
                fix(n);
                splayRake(rp.p);
            } else {
                int nch = N(n).ch[dir];
                pushNode(nch);
                N(n).rake = nch;  // demoted child becomes n's fresh rake root.
                N(nch).par = n;
                N(n).ch[dir] = node;
                N(node).par = n;
                fix(nch);
                fix(node);
                fix(n);
            }
            if (N(node).kind == EDGE) node = n;  // an edge leaf can't climb further; ride n.
        }
        return node;
    }

    int expose(int v) { return exposeRaw(handle_[static_cast<std::size_t>(v)]); }

    // put v and u on the two boundaries of one root cluster, v on side 0, u on side 1.
    // expose v, fence its root so exposing u can't tear it apart, then orient.
    void softExpose(int v, int u) {
        int root = expose(v);
        if (handle_[static_cast<std::size_t>(v)] == handle_[static_cast<std::size_t>(u)]) {
            if (endp(root, 1) == v || endp(root, 0) == u) {
                reverseNode(root);
                pushNode(root);
            }
            return;
        }
        if (N(root).kind == COMPRESS) {
            N(root).guard = true;
            int soot = expose(u);
            N(root).guard = false;
            pushNode(soot);
            fix(root);
            PD pd = parentDirComp(soot);
            if (pd.ok && pd.dir == 0) {
                reverseNode(root);
                pushNode(root);
            }
        }
    }

    // ---- link --------------------------------------------------------------
    // splice edge {v, u, weight}. mirror the reference exactly: expose each endpoint,
    // orient it, and hook the new edge either as a fresh compress child (endpoint is a
    // path end) or into a rake tree (endpoint is interior, so the branch is raked on).
    void linkRaw(int v, int u, F weight) {
        int hv = handle_[static_cast<std::size_t>(v)];
        int hu = handle_[static_cast<std::size_t>(u)];
        if (hv == -1 && hu == -1) {
            makeEdge(v, u, weight);
            return;
        }
        int e = makeEdge(v, u, weight);

        int left;  // the edge, married to u's side.
        if (hu == -1) {
            left = e;
        } else {
            int uu = exposeRaw(hu);
            pushNode(uu);
            if (endp(uu, 1) == u) {
                reverseNode(uu);
                pushNode(uu);
            }
            if (endp(uu, 0) == u) {  // u is a path end -- just compress the edge on.
                int nu = makeCompress(e, uu);
                N(e).par = nu;
                fix(e);
                N(uu).par = nu;
                fix(uu);
                fix(nu);
                left = nu;
            } else {  // u is interior -- its left branch gets raked onto u.
                int nu = uu;
                int leftCh = N(nu).ch[0];
                pushNode(leftCh);
                N(nu).ch[0] = e;
                N(e).par = nu;
                fix(e);
                int beta = N(nu).rake;
                int rk;
                if (beta != -1) {
                    pushNode(beta);
                    rk = makeRake(beta, leftCh);
                    N(beta).par = rk;
                    N(leftCh).par = rk;
                    fix(beta);
                    fix(leftCh);
                } else {
                    rk = leftCh;
                }
                fix(rk);
                N(nu).rake = rk;
                N(rk).par = nu;
                fix(rk);
                fix(nu);
                left = nu;
            }
        }

        if (hv == -1) return;  // v was fresh; `left` already pinned handle[v] via fix.
        int vv = exposeRaw(hv);
        pushNode(vv);
        if (endp(vv, 0) == v) {
            reverseNode(vv);
            pushNode(vv);
        }
        if (endp(vv, 1) == v) {  // v is a path end -- compress `left` on the right.
            int top = makeCompress(vv, left);
            N(vv).par = top;
            fix(vv);
            N(left).par = top;
            fix(left);
            fix(top);
        } else {  // v is interior -- its right branch gets raked onto v.
            int nv = vv;
            int rightCh = N(nv).ch[1];
            reverseNode(rightCh);
            pushNode(rightCh);
            N(nv).ch[1] = left;
            N(left).par = nv;
            fix(left);
            int alpha = N(nv).rake;
            int rk;
            if (alpha != -1) {
                pushNode(alpha);
                rk = makeRake(alpha, rightCh);
                N(alpha).par = rk;
                N(rightCh).par = rk;
                fix(alpha);
                fix(rightCh);
                fix(rk);
            } else {
                rk = rightCh;
            }
            N(nv).rake = rk;
            N(rk).par = nv;
            fix(rk);
            fix(nv);
        }
    }

    // ---- cut ---------------------------------------------------------------
    // root lost a child on one side; rebuild it from its rake tree -- promote the
    // rightmost raked subtree back to a path child, or, if the rake tree is empty,
    // dissolve the node and let the surviving child stand as a root.
    void bring(int root) {
        int rk = N(root).rake;
        if (rk == -1) {
            int left = N(root).ch[0];
            N(left).par = -1;
            fix(left);
        } else if (N(rk).kind != RAKE) {  // rake root is a single comp leaf.
            reverseNode(rk);
            pushNode(rk);
            N(root).ch[1] = rk;
            N(rk).par = root;
            N(root).rake = -1;
            fix(rk);
            fix(root);
        } else {  // walk to the rightmost rake leaf and splay it out.
            pushNode(rk);
            while (N(N(rk).ch[1]).kind == RAKE) {
                int right = N(rk).ch[1];
                pushNode(right);
                rk = right;
            }
            N(root).guard = true;
            splayRake(rk);
            N(root).guard = false;
            int newRake = N(rk).ch[0];
            int newRight = N(rk).ch[1];  // a comp leaf by construction.
            reverseNode(newRight);
            pushNode(newRight);
            N(root).ch[1] = newRight;
            N(newRight).par = root;
            N(root).rake = newRake;
            N(newRake).par = root;
            fix(newRake);
            fix(newRight);
            fix(root);
        }
    }

    void cutRaw(int v, int u) {
        softExpose(v, u);
        int root = handle_[static_cast<std::size_t>(v)];
        pushNode(root);
        // root is a compress with boundaries v..u; its right subtree carries u and the
        // edge {v,u}. detach it, flip it so the edge is exposed, and drop the edge from
        // both halves with bring().
        int right = N(root).ch[1];
        N(right).par = -1;
        reverseNode(right);
        pushNode(right);
        bring(right);
        bring(root);
    }

    // fold of the path v..u. after softExpose the root spans v..u; peel to the exact
    // path child depending on which boundaries the root actually exposes.
    F pathAggRaw(int v, int u) {
        softExpose(v, u);
        int root = handle_[static_cast<std::size_t>(v)];
        pushNode(root);
        if (endp(root, 0) == v && endp(root, 1) == u) return N(root).fold;
        if (endp(root, 0) == v) {
            pushNode(root);
            return N(N(root).ch[0]).fold;
        }
        if (endp(root, 1) == u) {
            pushNode(root);
            return N(N(root).ch[1]).fold;
        }
        pushNode(root);
        int n2 = N(root).ch[1];
        pushNode(n2);
        return N(N(n2).ch[0]).fold;
    }
};

// ===========================================================================
// the diameter cluster. len/d1/d2/diam exactly as documented at the top; a base
// edge of weight w seeds all four to w. compress and rake are the two merges from
// the header. reverse just swaps d1/d2 -- diam and len are reversal-symmetric.
// ===========================================================================
struct Diam {
    struct F {
        std::int64_t diam = 0;  // longest path fully inside the cluster.
        std::int64_t d1 = 0;    // max distance from boundary 1 to any inside vertex.
        std::int64_t d2 = 0;    // max distance from boundary 2 to any inside vertex.
        std::int64_t len = 0;   // weighted length of the boundary-to-boundary path.
    };
    static F id() { return F{0, 0, 0, 0}; }
    static F base(std::int64_t w) { return F{w, w, w, w}; }

    static F comp(const F& a, const F& b) {
        F r;
        r.diam = std::max({a.diam, b.diam, a.d2 + b.d1});
        r.d1 = std::max(a.d1, a.len + b.d1);
        r.d2 = std::max(b.d2, b.len + a.d2);
        r.len = a.len + b.len;
        return r;
    }
    // a is the retained path cluster; b hangs off a's boundary 2 (oriented so its own
    // boundary-2 faces the join). the header's rake formula, in code.
    static F rake(const F& a, const F& b) {
        F r;
        r.diam = std::max({a.diam, b.diam, a.d2 + b.d2});
        r.d1 = std::max(a.d1, a.len + b.d2);
        r.d2 = std::max(a.d2, b.d2);
        r.len = a.len;
        return r;
    }
    static void rev(F& x) { std::swap(x.d1, x.d2); }
};

// the problem's structure: link/cut a weighted forest, query a tree's diameter.
class DynamicDiameter {
public:
    explicit DynamicDiameter(int n) : f_(n) {}
    void link(int u, int v, std::int64_t w) { f_.link(u, v, w); }
    void cut(int u, int v) { f_.cut(u, v); }
    // diameter of the tree containing 1-based u -- 0 for a lone vertex.
    std::int64_t diameter(int u) { return f_.rootAgg(u).diam; }

private:
    TopForest<Diam> f_;
};

}  // namespace p0058
