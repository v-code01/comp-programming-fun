#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace p0028 {

// longest common *contiguous* substring of s and t.
//
// the whole idea in one line: glue the two strings with a separator that lives
// in neither -- w = s + '#' + t -- sort every suffix of w, and the answer falls
// out of two neighbors in that sorted order that came from different sides.
//
// why the separator. two suffixes share a prefix character-for-character until
// they diverge. an s-suffix hits '#' at position |s|; t never contains '#'. so
// no common prefix between an s-suffix and a t-suffix can reach across the
// boundary -- it's pinned inside s on one side and inside t on the other. that
// shared prefix *is* a common substring, sitting in both strings. drop the '#'
// and a suffix of s spills into t; the "match" would span the seam and lie about
// a substring that was never in s alone.
//
// why adjacent, different-half neighbors. in the suffix array the lcp of any two
// suffixes SA[i] and SA[j] (i<j) is the min of the lcp values strictly between
// them -- add a suffix in the middle and the shared prefix can only shrink. the
// best common substring is some cross pair, one suffix from each half. walk the
// sorted order from that pair's low end to its high end: the two ends live on
// different sides, so somewhere along the walk two *adjacent* suffixes straddle
// the seam. that adjacent lcp is at least the min over the whole span -- at least
// the pair's own lcp. so scanning adjacent cross pairs never misses the max… and
// it's the only place we have to look.
//
// construction: O(n log n) with n = |s| + |t| + 1. suffix array by prefix
// doubling -- sort by first char, then by the first 2, 4, 8… characters, each
// round a radix (counting) sort in O(n), log n rounds. kasai builds the lcp
// array in O(n). the adjacent walk is O(n). so O(n log n) end to end.
//
// lower bound. you have to read both strings -- Omega(|s| + |t|). the suffix
// array is a sort of n items; in the comparison model that's Omega(n log n). the
// doubling construction meets it. (SA-IS gets O(n) with a heavier machine; not
// worth it at n = 2e5.)

// sort the suffixes of `in` and return their start positions, smallest suffix
// first. a sentinel smaller than every real byte is appended so the sort of
// cyclic shifts of (in + sentinel) equals the suffix array -- the sentinel
// suffix lands at index 0 and gets dropped on the way out.
inline std::vector<int> suffix_array(const std::string& in) {
    std::string s = in;
    s.push_back('\1');  // sentinel -- below '#' (35) and every lowercase letter.
    const int n = static_cast<int>(s.size());
    const int alphabet = 256;  // one byte, treated unsigned.

    std::vector<int> p(n), c(n), cnt(std::max(alphabet, n), 0);

    // round zero: bucket-sort the single characters. `cnt` doubles as the class
    // table below, so it's sized for whichever is larger, the alphabet or n.
    for (int i = 0; i < n; ++i) ++cnt[static_cast<unsigned char>(s[i])];
    for (int i = 1; i < alphabet; ++i) cnt[i] += cnt[i - 1];
    for (int i = 0; i < n; ++i) p[--cnt[static_cast<unsigned char>(s[i])]] = i;
    c[p[0]] = 0;
    int classes = 1;
    for (int i = 1; i < n; ++i) {
        if (static_cast<unsigned char>(s[p[i]]) !=
            static_cast<unsigned char>(s[p[i - 1]]))
            ++classes;
        c[p[i]] = classes - 1;
    }

    std::vector<int> pn(n), cn(n);
    // each round sorts by the first 2^(h+1) characters, reusing last round's
    // classes as the key. the second half of the key is the class of the suffix
    // 2^h to the right -- a cyclic offset, which is why the sentinel matters.
    for (int h = 0; (1 << h) < n; ++h) {
        const int shift = 1 << h;
        // order by the second key first (a stable radix needs the low digit
        // sorted before the high one). pn holds that pre-order.
        for (int i = 0; i < n; ++i) {
            pn[i] = p[i] - shift;
            if (pn[i] < 0) pn[i] += n;
        }
        std::fill(cnt.begin(), cnt.begin() + classes, 0);
        for (int i = 0; i < n; ++i) ++cnt[c[pn[i]]];
        for (int i = 1; i < classes; ++i) cnt[i] += cnt[i - 1];
        // walk pn back-to-front so equal first-keys keep their second-key order.
        for (int i = n - 1; i >= 0; --i) p[--cnt[c[pn[i]]]] = pn[i];

        cn[p[0]] = 0;
        classes = 1;
        for (int i = 1; i < n; ++i) {
            const int cur1 = c[p[i]], cur2 = c[(p[i] + shift) % n];
            const int prv1 = c[p[i - 1]], prv2 = c[(p[i - 1] + shift) % n];
            if (cur1 != prv1 || cur2 != prv2) ++classes;
            cn[p[i]] = classes - 1;
        }
        c.swap(cn);
        if (classes == n) break;  // every suffix already unique -- done early.
    }

    // p[0] is the sentinel-only suffix, smaller than all. drop it and the rest
    // are the real suffixes of `in`, at their original positions.
    return std::vector<int>(p.begin() + 1, p.end());
}

// kasai's lcp array. lcp[i] is the longest common prefix of the two adjacent
// suffixes SA[i] and SA[i-1]; lcp[0] is 0 by convention. O(n): the match length
// h drops by at most one each step, so it never rewinds more than it advances.
inline std::vector<int> lcp_kasai(const std::string& s,
                                  const std::vector<int>& sa) {
    const int n = static_cast<int>(s.size());
    std::vector<int> rank(n, 0);
    for (int i = 0; i < n; ++i) rank[sa[i]] = i;

    std::vector<int> lcp(n, 0);
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (rank[i] > 0) {
            const int j = sa[rank[i] - 1];  // the suffix just before i's in SA.
            while (i + h < n && j + h < n && s[i + h] == s[j + h]) ++h;
            lcp[rank[i]] = h;
            if (h > 0) --h;  // next suffix loses its first char -- carry the rest.
        } else {
            h = 0;  // smallest suffix has no left neighbor to match against.
        }
    }
    return lcp;
}

// the answer. length of the longest substring common to s and t, 0 if none.
inline int longest_common_substring(const std::string& s, const std::string& t) {
    if (s.empty() || t.empty()) return 0;

    const int ns = static_cast<int>(s.size());
    std::string w;
    w.reserve(s.size() + t.size() + 1);
    w += s;
    w.push_back('#');  // separator -- not a lowercase letter, so no real match
    w += t;            // between the two halves can cross it.

    const std::vector<int> sa = suffix_array(w);
    const std::vector<int> lcp = lcp_kasai(w, sa);

    // which side of the seam a suffix starts on. position ns is the '#' itself --
    // it belongs to no real half, so it never forms a counted pair.
    auto half = [ns](int pos) -> int {
        if (pos < ns) return 0;   // s
        if (pos > ns) return 1;   // t
        return -1;                // the separator
    };

    int best = 0;
    for (int i = 1; i < static_cast<int>(sa.size()); ++i) {
        const int a = half(sa[i]);
        const int b = half(sa[i - 1]);
        // adjacent in sorted order, one from each half -- their shared prefix is
        // a real common substring. the longest such is the answer.
        if (a != -1 && b != -1 && a != b && lcp[i] > best) best = lcp[i];
    }
    return best;
}

}  // namespace p0028
