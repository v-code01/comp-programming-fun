#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace p0015 {

// codeforces 1009F -- dominant indices. rooted tree at vertex 1, n vertices.
// d(v, x) = how many vertices sit exactly x levels below v inside v's subtree.
// for each v report the x that maximizes d(v, x); ties break to the smallest x.
//
// the honest lower bound: you must read the tree (Omega(n)) and emit n answers
// (Omega(n)). so Omega(n) is the floor. this hits O(n) -- optimal to a constant.
//
// ---- the idea: long-path decomposition ----
//
// the depth histogram of v is d(v, .) with d(v, 0) = 1 (v itself) and
//        d(v, x) = sum over children c of d(c, x - 1).
// stack the children's histograms shifted right by one and add v's own 1 at the
// front. the naive merge is O(subtree size) per vertex -- O(n^2), quadratic on a
// path. the whole game is not paying for the merge.
//
// pick, at every vertex, its LONG child -- the child whose subtree reaches the
// most levels down. the long child's histogram, shifted by one, is almost all of
// the parent's histogram. so don't copy it: alias it. lay every histogram in one
// flat pool and give the long child the parent's block offset by +1. now
//        pool[off[v] + x]  ==  pool[off[son] + (x - 1)]      for x >= 1,
// for free -- same memory, no work. the parent writes its own 1 at off[v] and
// the long child's whole subtree is already accounted for one level down.
//
// only the SHORT children get merged by hand, level by level. that is the only
// arithmetic in the whole thing.
//
// ---- why the short merges total O(n) (the charging) ----
//
// follow long-child pointers from any vertex and you walk a "long path" down the
// tree. every vertex lies on exactly one long path. a short child c is the head
// of its own long path, and the cost of merging c is len[c] = the number of
// levels in c's subtree = the length of the long path that starts at c. sum that
// over all heads and you get the total length of all long paths -- which is n,
// since the paths partition the vertices. so the merges cost O(n). done.
//
// the same argument sizes the pool: fresh blocks are handed out only to long-path
// heads, len[head] cells each, summing to exactly n. one array of n ints.
//
// everything below is iterative. a path is a chain of depth 1e6 -- recursion
// would blow the stack -- so we work over an explicit BFS order, never the call
// stack.
//
// inputs: n >= 1; edges lists the n-1 undirected tree edges as 1-indexed pairs.
// output: vector of length n, out[i] = the dominant index of vertex i+1.
// time O(n), space O(n).
inline std::vector<int> dominant_indices(
    int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<int> ans(static_cast<std::size_t>(n < 0 ? 0 : n), 0);
    if (n <= 1) return ans;  // one vertex: only depth 0 exists, answer 0.

    const std::size_t N = static_cast<std::size_t>(n);

    // ---- CSR adjacency, undirected ----
    // count degrees, prefix-sum into head offsets, then scatter. one contiguous
    // neighbor array instead of n little vectors -- cache-friendly and a single
    // allocation. the whole build is O(n).
    std::vector<int> deg(N, 0);
    for (const auto& e : edges) {
        ++deg[static_cast<std::size_t>(e.first - 1)];
        ++deg[static_cast<std::size_t>(e.second - 1)];
    }
    std::vector<int> head(N + 1, 0);
    for (std::size_t v = 0; v < N; ++v) head[v + 1] = head[v] + deg[v];
    std::vector<int> adj(static_cast<std::size_t>(head[N]));
    std::vector<int> fill = head;  // mutable cursor per vertex.
    for (const auto& e : edges) {
        const int a = e.first - 1, b = e.second - 1;
        adj[static_cast<std::size_t>(fill[static_cast<std::size_t>(a)]++)] = b;
        adj[static_cast<std::size_t>(fill[static_cast<std::size_t>(b)]++)] = a;
    }

    // ---- BFS from root 0: parent, and an order with parent before child ----
    std::vector<int> parent(N, -1);
    std::vector<int> order;
    order.reserve(N);
    std::vector<char> seen(N, 0);
    order.push_back(0);
    seen[0] = 1;
    for (std::size_t qi = 0; qi < order.size(); ++qi) {
        const int u = order[qi];
        for (int i = head[static_cast<std::size_t>(u)];
             i < head[static_cast<std::size_t>(u) + 1]; ++i) {
            const int w = adj[static_cast<std::size_t>(i)];
            if (!seen[static_cast<std::size_t>(w)]) {
                seen[static_cast<std::size_t>(w)] = 1;
                parent[static_cast<std::size_t>(w)] = u;
                order.push_back(w);
            }
        }
    }

    // ---- len[] and son[], children before parents (reverse BFS) ----
    // len[u] = number of depth levels in u's subtree (leaf = 1). son[u] = the
    // long child, -1 if u is a leaf. reverse BFS visits every child before its
    // parent, so when we reach u its len is final and we push it up to parent.
    std::vector<int> len(N, 1);
    std::vector<int> son(N, -1);
    for (std::size_t k = order.size(); k-- > 0;) {
        const int u = order[k];
        const int p = parent[static_cast<std::size_t>(u)];
        // strict '>' keeps whichever long child we met first -- any longest works,
        // the tie among equal-height children never changes the answer.
        if (p >= 0 && len[static_cast<std::size_t>(u)] + 1 >
                          len[static_cast<std::size_t>(p)]) {
            len[static_cast<std::size_t>(p)] = len[static_cast<std::size_t>(u)] + 1;
            son[static_cast<std::size_t>(p)] = u;
        }
    }

    // ---- pool offsets, parents before children (forward BFS) ----
    // the long child aliases the parent's block at +1. every short child is a
    // long-path head and grabs a fresh len[c] block. the cursor lands exactly at
    // n -- that identity is the O(n) memory claim, made concrete.
    std::vector<int> off(N, 0);
    off[0] = 0;
    long long cursor = len[0];  // root is a head; its block is [0, len[0]).
    for (const int u : order) {
        const std::size_t su = static_cast<std::size_t>(u);
        if (son[su] >= 0) off[static_cast<std::size_t>(son[su])] = off[su] + 1;
        for (int i = head[su]; i < head[su + 1]; ++i) {
            const int w = adj[static_cast<std::size_t>(i)];
            if (w == parent[su] || w == son[su]) continue;  // up-edge or long child.
            off[static_cast<std::size_t>(w)] = static_cast<int>(cursor);
            cursor += len[static_cast<std::size_t>(w)];
        }
    }

    // ---- fill counts + compute answers, children before parents ----
    // counts stay <= n <= 1e6, so int holds every cell with room to spare.
    std::vector<int> pool(static_cast<std::size_t>(cursor), 0);
    for (std::size_t k = order.size(); k-- > 0;) {
        const int u = order[k];
        const std::size_t su = static_cast<std::size_t>(u);
        const int ou = off[su];
        pool[static_cast<std::size_t>(ou)] = 1;  // u sits alone at depth 0.

        if (son[su] < 0) {
            ans[su] = 0;  // leaf: depth 0 is the only depth, count 1.
            continue;
        }

        // the long child's histogram already lives at pool[ou+1 ..] via the alias
        // and was filled when the child was processed (it is deeper, so earlier in
        // this reverse pass). inherit its best depth, shifted down by one.
        int best_x = ans[static_cast<std::size_t>(son[su])] + 1;
        // if that winning count is only 1, then depth 0 -- also count 1 -- ties it
        // and wins by being smaller.
        if (pool[static_cast<std::size_t>(ou + best_x)] <= 1) best_x = 0;
        ans[su] = best_x;

        // merge the short children by hand. this is the only real work, and the
        // long-path charging caps the grand total at O(n).
        for (int i = head[su]; i < head[su + 1]; ++i) {
            const int w = adj[static_cast<std::size_t>(i)];
            if (w == parent[su] || w == son[su]) continue;
            const int ow = off[static_cast<std::size_t>(w)];
            const int lw = len[static_cast<std::size_t>(w)];  // lw <= len[u]-1.
            for (int j = 1; j <= lw; ++j) {
                // child depth j-1 lands at parent depth j.
                const int val = (pool[static_cast<std::size_t>(ou + j)] +=
                                 pool[static_cast<std::size_t>(ow + j - 1)]);
                const int cur_best = pool[static_cast<std::size_t>(ou + ans[su])];
                // strictly taller wins; equal height wins only if strictly nearer.
                if (val > cur_best || (val == cur_best && j < ans[su])) ans[su] = j;
            }
        }
    }

    return ans;
}

}  // namespace p0015
