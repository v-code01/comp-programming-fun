#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0036 {

// Dominator tree via Lengauer-Tarjan.
//
// w DOMINATES v if every path root->v runs through w. idom(v) is the deepest
// such w other than v -- the dominators of v form a chain from the root, and
// idom is the last link before v. we want idom for every reachable vertex.
//
// the naive test removes one vertex at a time and re-checks reachability --
// O(n*(m+n)). Lengauer-Tarjan does it in O(m*alpha(n)) with a path-compressed
// link-eval forest over the DFS tree. reading the graph is already Omega(n+m),
// so this is near-linear -- you can't beat the input by more than an inverse
// Ackermann.
//
// the whole trick is the SEMIDOMINATOR. number vertices by DFS preorder. define
//   sdom(w) = the vertex u of smallest dfs-number for which there is a path
//             u = v_0, v_1, ..., v_k = w  where every interior v_i (0<i<k) has
//             dfs-number GREATER than w's.
// intuition: sdom(w) is the highest point you can "jump" to w from while staying
// below w in the tree -- a lower bound on idom(w), and the pivot the two-pass
// idom derivation swings on. sdom is computed from predecessors:
//   sdom(w) = min over edges v->w of
//               dfn(v)                    if dfn(v) < dfn(w)   (a real jump up)
//               sdom(EVAL(v))             if dfn(v) > dfn(w)   (relayed down-tree)
// EVAL(v) returns, among the already-linked strict ancestors of v in the DFS
// forest, the one whose sdom has the smallest dfs-number -- that's the link-eval
// query, and path compression is what keeps it near-linear.
//
// idom falls out of sdom in two passes (Lengauer-Tarjan Theorem):
//   process w in REVERSE dfs order, bucket w under sdom(w). when the parent p of
//   w is linked in, drain bucket[p]: for each v there, let u = EVAL(v);
//   idom(v) = (sdom(u) < sdom(v)) ? u : p  -- either the immediate answer, or a
//   deferred pointer to be fixed up.
//   then process w in FORWARD dfs order: if idom(w) != sdom(w), idom(w) =
//   idom(idom(w)). one relaxation, ancestors already final.
//
// everything below lives in dfs-number space [1..cnt]; ord[] maps back to the
// caller's 1..n vertex ids at the end. root's idom is 0, and every unreachable
// vertex is 0.
inline std::vector<int> dominator_tree(
    int n, const std::vector<std::pair<int, int>>& edges, int root = 1) {
    std::vector<int> idom_out(static_cast<std::size_t>(n) + 1, 0);
    if (n <= 0) return idom_out;

    // ---- successor CSR. self-loops carry no reachability, so drop them here
    // once and never think about them again. ----
    const std::size_t m = edges.size();
    std::vector<int> sdeg(static_cast<std::size_t>(n) + 2, 0);
    for (const auto& e : edges) {
        if (e.first == e.second) continue;  // self-loop -- ignore.
        if (e.first < 1 || e.first > n || e.second < 1 || e.second > n) continue;
        ++sdeg[static_cast<std::size_t>(e.first) + 1];
    }
    for (int i = 1; i <= n + 1; ++i)
        sdeg[static_cast<std::size_t>(i)] +=
            sdeg[static_cast<std::size_t>(i) - 1];
    std::vector<int> shead = sdeg;  // shead[u]..shead[u+1] is u's out-edges.
    std::vector<int> sadj(m);
    for (const auto& e : edges) {
        if (e.first == e.second) continue;
        if (e.first < 1 || e.first > n || e.second < 1 || e.second > n) continue;
        sadj[static_cast<std::size_t>(shead[static_cast<std::size_t>(e.first)]++)] =
            e.second;
    }
    // shead got consumed as a cursor; rebuild the row starts from sdeg.
    shead = sdeg;

    // ---- iterative preorder DFS from root. a 2e5 chain would blow a recursive
    // stack, so walk it with an explicit one. dfn is 1-based; 0 means unseen. ----
    std::vector<int> dfn(static_cast<std::size_t>(n) + 1, 0);
    std::vector<int> ord(static_cast<std::size_t>(n) + 1, 0);   // dfn -> vertex
    std::vector<int> par(static_cast<std::size_t>(n) + 1, 0);   // dfn -> parent dfn
    int cnt = 0;

    std::vector<int> stk_v;   // current vertex on the walk
    std::vector<int> stk_it;  // its cursor into sadj
    stk_v.reserve(static_cast<std::size_t>(n));
    stk_it.reserve(static_cast<std::size_t>(n));

    dfn[static_cast<std::size_t>(root)] = ++cnt;
    ord[static_cast<std::size_t>(cnt)] = root;
    par[static_cast<std::size_t>(cnt)] = 0;
    stk_v.push_back(root);
    stk_it.push_back(shead[static_cast<std::size_t>(root)]);

    while (!stk_v.empty()) {
        const int u = stk_v.back();
        int& it = stk_it.back();
        const int stop = shead[static_cast<std::size_t>(u) + 1];
        bool descended = false;
        while (it < stop) {
            const int w = sadj[static_cast<std::size_t>(it++)];
            if (dfn[static_cast<std::size_t>(w)] == 0) {
                dfn[static_cast<std::size_t>(w)] = ++cnt;
                ord[static_cast<std::size_t>(cnt)] = w;
                par[static_cast<std::size_t>(cnt)] = dfn[static_cast<std::size_t>(u)];
                stk_v.push_back(w);
                stk_it.push_back(shead[static_cast<std::size_t>(w)]);
                descended = true;
                break;
            }
        }
        if (!descended) {
            stk_v.pop_back();
            stk_it.pop_back();
        }
    }

    if (cnt == 1) {
        idom_out[static_cast<std::size_t>(root)] = 0;  // lone reachable vertex.
        return idom_out;
    }

    // ---- predecessor CSR in dfs-number space. only reachable endpoints matter:
    // an unreachable source sits on no root-path and can't feed a semidominator.
    // ----
    std::vector<int> pdeg(static_cast<std::size_t>(cnt) + 2, 0);
    for (const auto& e : edges) {
        if (e.first == e.second) continue;
        const int du = (e.first >= 1 && e.first <= n)
                           ? dfn[static_cast<std::size_t>(e.first)]
                           : 0;
        const int dv = (e.second >= 1 && e.second <= n)
                           ? dfn[static_cast<std::size_t>(e.second)]
                           : 0;
        if (du == 0 || dv == 0) continue;
        ++pdeg[static_cast<std::size_t>(dv) + 1];
    }
    for (int i = 1; i <= cnt + 1; ++i)
        pdeg[static_cast<std::size_t>(i)] +=
            pdeg[static_cast<std::size_t>(i) - 1];
    std::vector<int> phead = pdeg;
    std::size_t pcount = static_cast<std::size_t>(pdeg[static_cast<std::size_t>(cnt) + 1]);
    std::vector<int> padj(pcount);
    for (const auto& e : edges) {
        if (e.first == e.second) continue;
        const int du = (e.first >= 1 && e.first <= n)
                           ? dfn[static_cast<std::size_t>(e.first)]
                           : 0;
        const int dv = (e.second >= 1 && e.second <= n)
                           ? dfn[static_cast<std::size_t>(e.second)]
                           : 0;
        if (du == 0 || dv == 0) continue;
        padj[static_cast<std::size_t>(phead[static_cast<std::size_t>(dv)]++)] = du;
    }
    phead = pdeg;

    // ---- link-eval state, all indexed by dfs-number. ----
    std::vector<int> semi(static_cast<std::size_t>(cnt) + 1);
    std::vector<int> label(static_cast<std::size_t>(cnt) + 1);
    std::vector<int> anc(static_cast<std::size_t>(cnt) + 1, 0);  // 0 == unlinked
    std::vector<int> idom(static_cast<std::size_t>(cnt) + 1, 0);
    for (int i = 1; i <= cnt; ++i) {
        semi[static_cast<std::size_t>(i)] = i;   // every vertex starts as its own sdom.
        label[static_cast<std::size_t>(i)] = i;
    }
    std::vector<std::vector<int>> bucket(static_cast<std::size_t>(cnt) + 1);

    // iterative path compression -- a long ancestor chain must not recurse. we
    // collect the nodes whose grandparent is still linked, deepest last, then
    // fold labels from the deepest one down so each sees a finalized parent.
    std::vector<int> comp;  // scratch reused across calls.
    auto eval = [&](int v) -> int {
        if (anc[static_cast<std::size_t>(v)] == 0)
            return label[static_cast<std::size_t>(v)];
        comp.clear();
        int u = v;
        while (anc[static_cast<std::size_t>(anc[static_cast<std::size_t>(u)])] != 0) {
            comp.push_back(u);
            u = anc[static_cast<std::size_t>(u)];
        }
        for (std::size_t k = comp.size(); k-- > 0;) {
            const int x = comp[k];
            const int a = anc[static_cast<std::size_t>(x)];
            if (semi[static_cast<std::size_t>(label[static_cast<std::size_t>(a)])] <
                semi[static_cast<std::size_t>(label[static_cast<std::size_t>(x)])])
                label[static_cast<std::size_t>(x)] = label[static_cast<std::size_t>(a)];
            anc[static_cast<std::size_t>(x)] = anc[static_cast<std::size_t>(a)];
        }
        return label[static_cast<std::size_t>(v)];
    };

    // ---- pass 1: reverse dfs order. compute sdom, bucket, link, drain. ----
    for (int w = cnt; w >= 2; --w) {
        const std::size_t wsz = static_cast<std::size_t>(w);
        for (int e = phead[static_cast<std::size_t>(w)];
             e < phead[static_cast<std::size_t>(w) + 1]; ++e) {
            const int v = padj[static_cast<std::size_t>(e)];
            const int u = eval(v);
            if (semi[static_cast<std::size_t>(u)] < semi[wsz])
                semi[wsz] = semi[static_cast<std::size_t>(u)];
        }
        bucket[static_cast<std::size_t>(semi[wsz])].push_back(w);

        const int p = par[wsz];
        anc[wsz] = p;  // LINK(par[w], w): w now hangs off its tree parent.

        auto& pb = bucket[static_cast<std::size_t>(p)];
        for (const int v : pb) {
            const int u = eval(v);
            // v's sdom is p (that's why it sat in p's bucket). if the best
            // ancestor u has a shallower sdom, idom(v) relays to u; else it's p.
            idom[static_cast<std::size_t>(v)] =
                (semi[static_cast<std::size_t>(u)] <
                 semi[static_cast<std::size_t>(v)])
                    ? u
                    : p;
        }
        pb.clear();
    }

    // ---- pass 2: forward dfs order. resolve the relayed pointers. ----
    for (int w = 2; w <= cnt; ++w) {
        if (idom[static_cast<std::size_t>(w)] != semi[static_cast<std::size_t>(w)])
            idom[static_cast<std::size_t>(w)] =
                idom[static_cast<std::size_t>(idom[static_cast<std::size_t>(w)])];
    }
    idom[1] = 0;  // root has no dominator.

    // ---- back to the caller's vertex ids. unreachable vertices stay 0. ----
    for (int i = 2; i <= cnt; ++i)
        idom_out[static_cast<std::size_t>(ord[static_cast<std::size_t>(i)])] =
            ord[static_cast<std::size_t>(idom[static_cast<std::size_t>(i)])];
    idom_out[static_cast<std::size_t>(root)] = 0;
    return idom_out;
}

}  // namespace p0036
