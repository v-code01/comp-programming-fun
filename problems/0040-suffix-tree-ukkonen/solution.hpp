#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace p0040 {

// two numbers off one structure:
//   (1) how many distinct non-empty substrings s has,
//   (2) the length of the longest substring that shows up at least twice.
//
// s can hold ~2e10 distinct substrings at |s| = 2e5. we are never going to
// enumerate those. the suffix tree holds all of them in O(n) nodes -- every
// substring of s is a prefix of some suffix, every suffix is a root-to-leaf
// path, so every substring is a *position on an edge*. count the positions, not
// the strings. that compression is the whole trick, and it is why (1) is a sum
// of edge lengths instead of a set.
//
// the tree: a compressed trie of every suffix. edges carry a slice [start, end)
// of the text -- two ints, not a string -- so the whole thing is O(n) memory.
//
// sentinel. append a character that appears nowhere in s ('{' in effect, code
// 26 here). then no suffix is a prefix of another, so every suffix ends at its
// own leaf and every internal node is a real branch. without it a suffix like
// "a" in "aa" dies in the middle of an edge and the counting falls apart.
//
// ---- ukkonen, online, left to right ----
//
// after reading s[0..i) the tree is the suffix tree of that prefix. read s[i]
// and every suffix has to grow by one character. the three moves:
//
//   rule 1 -- a suffix that ends at a leaf. its edge just got one longer. don't
//     touch it. leaves store end = -1 and read a single global `leaf_end_`, so
//     one assignment per step extends all of them at once. O(1), not O(leaves).
//
//   rule 2 -- the suffix ends somewhere with no s[i] continuation. hang a new
//     leaf. if it ended mid-edge, split the edge first: a fresh internal node
//     at the split point, the old edge's tail below it, the new leaf beside.
//
//   rule 3 -- the suffix ends somewhere that *already* continues with s[i].
//     nothing to do… and the show stopper: if this suffix is already in the
//     tree, so is every shorter one (they're suffixes of it). stop the phase
//     right there. that early exit is what keeps the work linear.
//
// the active point (active_node, active_edge, active_len) is where the next
// suffix to be inserted currently ends: from active_node, take the edge whose
// first character is text[active_edge], walk active_len characters down it.
// `remainder_` counts the suffixes still owed an insert.
//
// suffix links. after handling suffix c*alpha (c a character, alpha a string),
// the next one to handle is alpha. an internal node for c*alpha links to the
// node for alpha, so we hop there in O(1) instead of re-walking from the root.
// the link of a node created by a split is filled in on the *next* split (or on
// the next rule-2/rule-3 event) inside the same phase -- that's what
// `last_internal` is carrying around.
//
// ---- complexity ----
//
// O(n * ALPHA) time, O(n * ALPHA) memory, ALPHA = 27 (children as a flat array,
// so a child lookup is one load). the amortization: active_len_ goes up by at
// most 1 per phase, and every walk-down or suffix-link hop only ever decreases
// it. it can't fall more than it rose, so the total hopping across the whole
// build telescopes to O(n). nodes: each phase creates at most one leaf per
// insert and each split adds one internal node, and there are exactly n+1
// leaves -- so < 2n+2 nodes.
//
// lower bound, honestly. you must read s: Omega(|s|). ukkonen matches it. and
// note what's *not* a lower bound here -- the number of distinct substrings is
// Theta(n^2) for a random string, so any answer to (1) that enumerates them is
// quadratic by construction. the tree beats that bound because it never
// materializes a substring: it counts edge positions. n nodes, n^2 substrings.
//
// ---- reading the answers off ----
//
// (1) every distinct non-empty substring of the text is exactly one position on
//     exactly one edge (walk the path, stop anywhere). so sum the edge lengths.
//     substrings *containing* the sentinel are not substrings of s -- and since
//     the sentinel is the last character of the text and occurs once, an edge
//     can only ever contain it as its final character. so drop one character
//     from every edge that runs to the end of the text. the sentinel-only edge
//     off the root contributes 0, exactly as it should.
//
// (2) an internal node is a branching node: two different characters follow it,
//     so its path occurs at at least two positions in the text. conversely any
//     substring that occurs twice is a prefix of some branching path. so the
//     longest repeated substring is the string-depth of the deepest internal
//     node. no internal node (every character distinct) -> 0. overlapping
//     occurrences count, which is what the tree gives you for free -- "ana" in
//     "banana" is two overlapping hits and one internal node.
//     the only internal nodes are the root and the splits, and a split's depth
//     is the active point's depth at the moment it's cut -- known, right then,
//     for free. write it down and answer (2) is a max over an int array. no
//     traversal at all: a dfs over 4e5 nodes costs a cache miss per hop and was
//     measurably as expensive as the entire build.
class SuffixTree {
public:
    // ALPHA-1 lowercase letters plus the sentinel in the last slot.
    static constexpr int kAlpha = 27;
    static constexpr std::uint8_t kSentinel = 26;

    explicit SuffixTree(const std::string& s) {
        text_.reserve(s.size() + 1);
        for (char c : s) {
            text_.push_back(static_cast<std::uint8_t>(c - 'a'));
        }
        text_.push_back(kSentinel);
        m_ = static_cast<int>(text_.size());

        // hard ceiling: n+1 leaves and at most n+1 internal nodes, plus root.
        // reserve once -- the build never reallocates, so raw indices stay put.
        nodes_.reserve(2 * static_cast<std::size_t>(m_) + 2);
        kids_.reserve(static_cast<std::size_t>(m_) + 2);
        idepth_.reserve(static_cast<std::size_t>(m_) + 2);
        make_internal(0, 0, 0);  // root. index 0. its link points at itself.

        for (int i = 0; i < m_; ++i) extend(i);
    }

    // {distinct substrings of s, length of the longest substring occurring >= 2}
    //
    // two flat scans over arrays. no tree walk -- a dfs over 4e5 nodes chases
    // pointers through tens of megabytes and eats a cache miss per hop; it cost
    // as much as the whole build. neither answer needs the shape of the tree,
    // only its edge lengths and its internal depths, and both of those are
    // already sitting in memory, in order.
    std::pair<std::int64_t, std::int64_t> stats() const {
        // (1) sum of edge lengths, sentinel character excluded. every node but
        // the root owns exactly one incoming edge, so the node array *is* the
        // edge list.
        std::int64_t distinct = 0;
        for (std::size_t v = 1; v < nodes_.size(); ++v) {
            distinct += label_len(static_cast<int>(v));
        }

        // (2) the deepest internal node. every internal node's string depth was
        // written down when the split created it -- and a node's depth never
        // moves after that, because its path label doesn't. so just take the
        // max. the root sits at 0, which is also the right answer when there is
        // no internal node at all -- every substring unique, nothing repeats.
        std::int64_t best = 0;
        for (int d : idepth_) best = std::max(best, static_cast<std::int64_t>(d));
        return {distinct, best};
    }

private:
    struct Node {
        int start;  // first text index of the incoming edge label.
        int end;    // one past the last. -1 means "the global leaf end".
        int link;   // suffix link. root for anything that never got one.
        int kids;   // block index into kids_. -1 for a leaf -- see below.
    };

    // leaves hold no child block. a leaf in ukkonen stays a leaf forever: a
    // split cuts the *edge above* a node and inserts a fresh internal node, it
    // never grows children under the old one. so half the nodes never need 27
    // slots, and not allocating them halves the tree's memory and keeps the hot
    // node array small enough to stay in cache.
    int child(int v, int c) const {
        return kids_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)]
                                                  .kids)]
                    [static_cast<std::size_t>(c)];
    }

    void set_child(int v, int c, int u) {
        kids_[static_cast<std::size_t>(nodes_[static_cast<std::size_t>(v)].kids)]
             [static_cast<std::size_t>(c)] = u;
    }

    int make_leaf(int start) {
        nodes_.push_back(Node{start, -1, 0, -1});
        return static_cast<int>(nodes_.size()) - 1;
    }

    // `depth` is the length of the node's path from the root -- its string
    // depth, recorded once, at O(1), because that's answer (2) and reading it
    // back off the tree later would cost a traversal.
    int make_internal(int start, int end, int depth) {
        std::array<int, kAlpha> empty;
        empty.fill(-1);
        kids_.push_back(empty);
        idepth_.push_back(depth);
        nodes_.push_back(Node{start, end, 0,
                              static_cast<int>(kids_.size()) - 1});
        return static_cast<int>(nodes_.size()) - 1;
    }

    int depth_of(int v) const {
        return idepth_[static_cast<std::size_t>(
            nodes_[static_cast<std::size_t>(v)].kids)];
    }

    // leaves don't store their end -- they read the global one. that's rule 1:
    // one write per phase grows every leaf.
    int edge_end(int v) const {
        const int e = nodes_[static_cast<std::size_t>(v)].end;
        return e < 0 ? leaf_end_ : e;
    }

    // edge length with the sentinel character not counted. the sentinel sits at
    // text_[m_-1] and nowhere else, so an edge touches it iff it runs to m_.
    std::int64_t label_len(int v) const {
        const int e = edge_end(v);
        std::int64_t len = e - nodes_[static_cast<std::size_t>(v)].start;
        if (e == m_) --len;
        return len;
    }

    void extend(int pos) {
        leaf_end_ = pos + 1;  // rule 1, for the whole leaf set, in one store.
        ++remainder_;         // one more suffix owed.
        int last_internal = -1;  // needs a suffix link before the phase ends.

        while (remainder_ > 0) {
            if (active_len_ == 0) active_edge_ = pos;  // sitting on the node.

            // active_node_ is the root or a split -- always internal, so it has
            // a child block. it can't be a leaf: a leaf edge runs to pos+1 while
            // the active string only reaches pos, so the walk-down below never
            // has the length to step into one.
            const int c = text_[static_cast<std::size_t>(active_edge_)];
            const int nxt = child(active_node_, c);

            if (nxt == -1) {
                // rule 2, no split -- nothing leaves active_node_ on c.
                set_child(active_node_, c, make_leaf(pos));
                if (last_internal != -1) {
                    nodes_[static_cast<std::size_t>(last_internal)].link =
                        active_node_;
                    last_internal = -1;
                }
            } else {
                const int len =
                    edge_end(nxt) - nodes_[static_cast<std::size_t>(nxt)].start;
                if (active_len_ >= len) {
                    // the active point ran past this edge. skip to the node
                    // below and re-read. this is the step that pays for itself:
                    // it only ever *shrinks* active_len_.
                    active_edge_ += len;
                    active_len_ -= len;
                    active_node_ = nxt;
                    continue;
                }

                const int at = nodes_[static_cast<std::size_t>(nxt)].start +
                               active_len_;
                if (text_[static_cast<std::size_t>(at)] ==
                    text_[static_cast<std::size_t>(pos)]) {
                    // rule 3. the character is already on the edge -- this
                    // suffix is present, and so is every shorter one. link up
                    // whatever split is still waiting, step one deeper, and
                    // stop the phase. show stopper.
                    if (last_internal != -1 && active_node_ != 0) {
                        nodes_[static_cast<std::size_t>(last_internal)].link =
                            active_node_;
                        last_internal = -1;
                    }
                    ++active_len_;
                    break;
                }

                // rule 2 with a split. cut `nxt` at active_len_: a new internal
                // node takes the head, `nxt` keeps the tail, and the new leaf
                // hangs off the side with text[pos].
                // the split sits active_len_ characters below active_node_ --
                // that's the active point, and that's its string depth.
                const int split = make_internal(
                    nodes_[static_cast<std::size_t>(nxt)].start, at,
                    depth_of(active_node_) + active_len_);
                set_child(active_node_, c, split);
                set_child(split, text_[static_cast<std::size_t>(pos)],
                          make_leaf(pos));

                nodes_[static_cast<std::size_t>(nxt)].start = at;
                set_child(split, text_[static_cast<std::size_t>(at)], nxt);

                // the previous split in this phase was the same string with one
                // more character in front. that's exactly a suffix link.
                if (last_internal != -1) {
                    nodes_[static_cast<std::size_t>(last_internal)].link = split;
                }
                last_internal = split;
            }

            --remainder_;  // that suffix is in. on to the next shorter one.

            if (active_node_ == 0 && active_len_ > 0) {
                // at the root there's no link to follow: drop the first
                // character by hand. remainder_ tells us which suffix is next.
                --active_len_;
                active_edge_ = pos - remainder_ + 1;
            } else if (active_node_ != 0) {
                // c*alpha handled, alpha is next -- one hop, no re-walk.
                active_node_ = nodes_[static_cast<std::size_t>(active_node_)].link;
            }
        }
    }

    std::vector<std::uint8_t> text_;  // s mapped to 0..25, then the sentinel.
    std::vector<Node> nodes_;                    // 16 bytes each, leaves included.
    std::vector<std::array<int, kAlpha>> kids_;  // one block per internal node.
    std::vector<int> idepth_;                    // string depth, same indexing.
    int m_ = 0;

    // the active point plus the debt.
    int active_node_ = 0;
    int active_edge_ = 0;
    int active_len_ = 0;
    int remainder_ = 0;
    int leaf_end_ = 0;
};

// build the tree, read both numbers off it. O(n * ALPHA).
inline std::pair<std::int64_t, std::int64_t> suffix_tree_stats(
    const std::string& s) {
    return SuffixTree(s).stats();
}

}  // namespace p0040
