#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace p0012 {

// count, for each query x, how many length-|x| windows of s (counted by
// position) equal *some* rotation of x. a position that matches several
// rotations is still one position -- so we count each matched position once.
//
// the whole trick: the set of rotations of x is exactly the set of length-|x|
// windows of the doubled string x+x. a position of s is a specific string; it
// equals at most one rotation-string (two distinct rotation-*strings* can't both
// equal it). so if we knew, for every distinct rotation-string, how many
// positions of s spell it, the answer is the sum -- with no cross-position
// double counting. the suffix automaton of s gives exactly that: each substring
// lives in one state, and that state's endpos size is the number of positions
// where it ends.
//
// build the SAM of s once, tag every state with its endpos size (occurrence
// count), then for each query walk x+x through the automaton keeping the current
// matched suffix at length |x|. every time the walk sits on a state spelling a
// length-|x| rotation, add that state's count -- but only the first time we land
// on that state this query, so a rotation that repeats (think "aa") or two
// offsets that spell the same string don't get counted twice.
//
// build: O(|s| * A) with A = 26. query: O(|x|) amortized -- the doubled walk is
// linear because suffix-link fallbacks are paid for by the length we advance.
// total O(|s| * A + sum|x|). space O(|s| * A) for the transition table.
//
// lower bound: you must read s and every x, so Omega(|s| + sum|x|). the SAM hits
// that up to the alphabet factor A in construction and the transition table.
class SuffixAutomaton {
public:
    static constexpr int kAlpha = 26;

    // build the automaton of s and precompute endpos sizes.
    void build(const std::string& s) {
        n_ = static_cast<int>(s.size());
        const int cap = 2 * n_ + 5;  // a SAM over n chars has < 2n states.
        len_.assign(cap, 0);
        link_.assign(cap, -1);
        cnt_.assign(cap, 0);
        next_.assign(cap, std::array<int, kAlpha>{});  // 0 == "no edge".
        sz_ = 0;
        last_ = new_state(0);  // the initial state -- id 0, empty string, len 0.
        // link_[0] stays -1. it is the only state with len 0, so it's also the
        // unique root when we bucket-sort by len below.
        for (char ch : s) extend(ch - 'a');
        compute_counts();
        visited_.assign(sz_, 0);  // per-query dedupe marks, 0 == "never counted".
    }

    // answer one query. query_id must be a nonzero value unique to this call
    // over the lifetime of this automaton -- it's the dedupe tag written into
    // visited_, and 0 means "unmarked", so ids start at 1.
    std::int64_t count_cyclic(const std::string& x, int query_id) {
        const int L = static_cast<int>(x.size());
        if (L == 0) return 0;

        int v = 0;   // current state.
        int l = 0;   // length of the longest suffix matched so far, capped at L.
        std::int64_t ans = 0;

        // rotation k of x is the length-L window of x+x starting at offset k,
        // k = 0..L-1. that window *ends* at index k+L-1, so ends range over
        // [L-1, 2L-2]. process 2L-1 characters of the doubled string and every
        // ending index is covered exactly once.
        const int steps = 2 * L - 1;
        for (int i = 0; i < steps; ++i) {
            const int c = x[i % L] - 'a';

            if (next_[v][c] != 0) {
                // the last (l+1) characters are still a substring of s -- extend.
                v = next_[v][c];
                ++l;
            } else {
                // dead end. drop to the longest suffix that *can* eat c by
                // walking suffix links. each link strictly shortens l, and l only
                // ever grew by +1 per character -- that's why the walk is linear.
                while (v != -1 && next_[v][c] == 0) v = link_[v];
                if (v == -1) {
                    v = 0;  // nothing matches c; restart from the empty suffix.
                    l = 0;
                } else {
                    l = len_[v] + 1;  // matched suffix is now (this state)+c.
                    v = next_[v][c];
                }
            }

            // cap the match at L. a rotation is exactly L long, and the length-L
            // suffix depends only on the last L characters -- so throwing away
            // anything longer costs nothing, and it keeps v on the state that
            // actually spells the length-L window. edge next_[v][c] is defined
            // iff (any string of v)+c is a substring, a property of v's endpos
            // set shared by every string in v -- so extending the capped state
            // still tracks the true length-L suffix.
            if (l > L) {
                while (len_[link_[v]] >= L) v = link_[v];
                l = L;
            }

            // sitting on a length-L rotation that occurs in s. add its positions
            // once: distinct rotation-strings live on distinct states, so marking
            // the state dedupes repeated rotations without dropping real ones.
            if (i >= L - 1 && l == L) {
                if (visited_[v] != query_id) {
                    visited_[v] = query_id;
                    ans += cnt_[v];
                }
            }
        }
        return ans;
    }

    int states() const { return sz_; }

private:
    int new_state(int length) {
        const int id = sz_++;
        len_[id] = length;
        link_[id] = -1;
        cnt_[id] = 0;  // clones keep this 0; real states get 1 in extend().
        next_[id].fill(0);
        return id;
    }

    // the online SAM extension (cp-algorithms / Blumer). append one character and
    // splice in at most one new state plus at most one clone.
    void extend(int c) {
        const int cur = new_state(len_[last_] + 1);
        cnt_[cur] = 1;  // a real prefix-ending state -- it owns one endpos.

        int p = last_;
        // add the c-edge to cur along the suffix-link chain until some state
        // already has a c-edge.
        while (p != -1 && next_[p][c] == 0) {
            next_[p][c] = cur;
            p = link_[p];
        }
        if (p == -1) {
            link_[cur] = 0;  // c never appeared before -- link to root.
        } else {
            const int q = next_[p][c];
            if (len_[p] + 1 == len_[q]) {
                link_[cur] = q;  // q is a clean continuation -- reuse it.
            } else {
                // q packs strings longer than len_[p]+1. split off a clone that
                // holds exactly the length-(len_[p]+1) equivalence class, and it
                // inherits q's edges but owns no new endpos -- cnt stays 0.
                const int clone = new_state(len_[p] + 1);
                next_[clone] = next_[q];
                link_[clone] = link_[q];
                while (p != -1 && next_[p][c] == q) {
                    next_[p][c] = clone;
                    p = link_[p];
                }
                link_[q] = link_[cur] = clone;
            }
        }
        last_ = cur;
    }

    // endpos size of a state = number of positions in s where its strings end.
    // real states start at 1 (their own prefix), clones at 0, then counts flow up
    // the suffix-link tree: a suffix's positions are the union of its extensions'.
    // process states in decreasing len so a child is always added before its
    // (shorter) parent. bucket sort by len -- len in [0, n], so it's O(n).
    void compute_counts() {
        std::vector<int> bucket(n_ + 1, 0);
        for (int i = 0; i < sz_; ++i) ++bucket[len_[i]];
        for (int i = 1; i <= n_; ++i) bucket[i] += bucket[i - 1];
        std::vector<int> order(sz_);
        for (int i = 0; i < sz_; ++i) order[--bucket[len_[i]]] = i;
        // order is ascending by len; order[0] is the root (the only len-0 state).
        // walk it backwards -- longest first -- and push each count to its link.
        for (int i = sz_ - 1; i >= 1; --i) {
            const int v = order[i];
            if (link_[v] >= 0) cnt_[link_[v]] += cnt_[v];
        }
    }

    int n_ = 0;
    int sz_ = 0;
    int last_ = 0;
    std::vector<int> len_;
    std::vector<int> link_;
    std::vector<int> cnt_;                        // endpos sizes.
    std::vector<std::array<int, kAlpha>> next_;   // 0 sentinel == no transition.
    std::vector<int> visited_;                    // last query_id that counted a state.
};

// convenience: build the automaton for s once and answer every query.
inline std::vector<std::int64_t> solve(const std::string& s,
                                       const std::vector<std::string>& queries) {
    SuffixAutomaton sam;
    sam.build(s);
    std::vector<std::int64_t> out;
    out.reserve(queries.size());
    for (std::size_t i = 0; i < queries.size(); ++i) {
        out.push_back(sam.count_cyclic(queries[i], static_cast<int>(i) + 1));
    }
    return out;
}

}  // namespace p0012
