#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace p0030 {

// count the distinct palindromic substrings of s. the empty string doesn't count.
//
// the tool is the palindromic tree -- the eertree. every node is one distinct
// palindrome. so the answer is (number of nodes) minus the two roots. build it in
// one left-to-right pass and read off the size.
//
// -- the two roots --
// the tree starts seeded with two imaginary nodes, and getting these right is the
// whole game.
//   node 0 -- length -1. the odd root. it exists so a single character has a
//     parent. a length-1 palindrome "extends" the length -1 root by wrapping one
//     char on each side: len goes -1 -> +1. the character-before test for it,
//     s[i - (-1) - 1] == s[i], reads s[i] == s[i] -- true for free, always. that's
//     why the walk can never fall off the bottom.
//   node 1 -- length 0. the even root. the empty string. every even-length
//     palindrome descends from here.
//   link(0) = 0, link(1) = 0. both point at the odd root.
// neither root is a real substring, so the final count subtracts 2.
//
// -- adding a character --
// we carry `last`, the node for the longest palindromic suffix of the prefix so
// far. to extend by s[i] we walk suffix links up from `last` until the palindrome
// at the node can be wrapped in s[i] -- that is, until s[i - len(v) - 1] == s[i].
// call that node t. the new longest palindromic suffix is t with an s[i] glued on
// each end, length len(t) + 2.
//   if the edge t --s[i]--> already exists, that palindrome is already in the
//   tree. nothing new. just move `last` onto it.
//   otherwise it's new: one fresh node. its suffix link is found by continuing the
//   same walk from link(t) to the next node that s[i] can wrap, then taking that
//   node's s[i] edge -- which is guaranteed to exist, it's a shorter suffix we
//   built earlier. the single length-1 node links straight to the even root.
//
// -- why it's O(|s|) amortized --
// each character adds at most one node: at most one new distinct palindromic
// suffix ends at each position. the cost of a step is the length of the suffix-link
// walk. that walk climbs, and the only thing that pushes `last` back down is
// adding a node -- so the total climb across the whole string telescopes, the same
// argument that makes manacher's linear. total link-walk work is O(|s|). the
// alphabet factor lives only in the per-node transition array, O(|s| * 26) space.
//
// -- lower bound, honest --
// reading s is Omega(|s|). and you can't beat linear on the output either: the
// count can be as large as |s| itself. "aaaa" gives a, aa, aaa, aaaa -- one
// distinct palindrome per length, |s| of them. so no algorithm reports the count
// in sub-linear time. the eertree is optimal.
class Eertree {
public:
    // seed the two roots. reserve for n more nodes -- at most one per character.
    explicit Eertree(std::size_t n) {
        len_.reserve(n + 2);
        link_.reserve(n + 2);
        next_.reserve((n + 2) * kAlpha);
        str_.reserve(n);
        new_node(-1);   // node 0 -- odd root, length -1.
        new_node(0);    // node 1 -- even root, length 0.
        link_[0] = 0;
        link_[1] = 0;
        last_ = 1;      // longest palindromic suffix of the empty prefix: nothing.
    }

    // append one letter c in [0, 26). grows the tree by at most one node.
    void add(int c) {
        str_.push_back(static_cast<char>('a' + c));
        const int i = static_cast<int>(str_.size()) - 1;

        // walk up to the longest palindromic suffix that s[i] can wrap.
        const int t = extendable(last_, i);
        if (next_[static_cast<std::size_t>(t) * kAlpha + c] != -1) {
            last_ = next_[static_cast<std::size_t>(t) * kAlpha + c];
            return;  // this palindrome already lives in the tree.
        }

        const int cur = new_node(len_[t] + 2);  // may reallocate -- read after.
        if (len_[cur] == 1) {
            link_[cur] = 1;  // a lone char hangs off the even root.
        } else {
            // keep climbing from t's link for the suffix link's target node.
            const int u = extendable(link_[t], i);
            link_[cur] = next_[static_cast<std::size_t>(u) * kAlpha + c];
        }
        next_[static_cast<std::size_t>(t) * kAlpha + c] = cur;
        last_ = cur;
    }

    // the two roots aren't real substrings, so the palindrome count drops them.
    std::int64_t distinct_count() const {
        return static_cast<std::int64_t>(len_.size()) - 2;
    }

private:
    static constexpr std::size_t kAlpha = 26;

    std::vector<int> len_;    // palindrome length per node.
    std::vector<int> link_;   // suffix link per node.
    std::vector<int> next_;   // flat transition table: node * 26 + c -> node/-1.
    std::string str_;         // characters appended so far.
    int last_ = 1;

    int new_node(int length) {
        len_.push_back(length);
        link_.push_back(0);
        next_.insert(next_.end(), kAlpha, -1);
        return static_cast<int>(len_.size()) - 1;
    }

    // climb suffix links from v until the palindrome there can be wrapped in the
    // char at position i -- s[i - len(v) - 1] == s[i]. the odd root's -1 length
    // makes that test s[i] == s[i], so this always terminates.
    int extendable(int v, int i) const {
        while (true) {
            const int j = i - len_[v] - 1;
            if (j >= 0 && str_[static_cast<std::size_t>(j)] ==
                              str_[static_cast<std::size_t>(i)]) {
                return v;
            }
            v = link_[v];
        }
    }
};

// distinct nonempty palindromic substrings of s, in O(|s| * 26).
inline std::int64_t count_distinct(const std::string& s) {
    Eertree tree(s.size());
    for (const char c : s) {
        tree.add(static_cast<unsigned char>(c) - 'a');
    }
    return tree.distinct_count();
}

}  // namespace p0030
