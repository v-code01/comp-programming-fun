#include <cstdint>
#include <cstdio>
#include <string>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// cpf::Rng deals in ints; the problem eats strings. map a base-alpha digit onto
// a lowercase alphabet, and keep the alphabet small on purpose -- 2 or 3 letters
// so common substrings actually happen and the answer isn't always 1.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    return s;
}

// one pair -- the differential's "solution" side.
struct Pair {
    std::string s;
    std::string t;
    bool operator==(const Pair&) const = default;
};

}  // namespace

int main() {
    using namespace p0028;

    // ---- hand examples, ground truth ----
    CPF_EQ(longest_common_substring("abcde", "cdef"), 3);   // "cde".
    CPF_EQ(longest_common_substring("abc", "xyz"), 0);      // disjoint alphabets.
    CPF_EQ(longest_common_substring("abcabc", "abcabc"), 6);  // identical.
    CPF_EQ(longest_common_substring("abcdef", "cd"), 2);   // t is a substring of s.
    CPF_EQ(longest_common_substring("cd", "abcdef"), 2);   // and the other way.
    CPF_EQ(longest_common_substring("a", "a"), 1);         // shortest match.
    CPF_EQ(longest_common_substring("a", "b"), 0);         // shortest miss.
    CPF_EQ(longest_common_substring("aaaaa", "aaa"), 3);   // run capped by |t|.
    CPF_EQ(longest_common_substring("banana", "ananas"), 5);  // "anana".
    CPF_EQ(longest_common_substring("xabcx", "yabcy"), 3);  // match in the middle.

    // the separator has to stop a match from crossing the seam. these two share
    // no letter, so any answer above 0 means s spilled into t through the '#'.
    CPF_EQ(longest_common_substring("abc", "def"), 0);
    // last char of s equals first char of t -- a single letter, not a bridge.
    CPF_EQ(longest_common_substring("abcx", "xdef"), 1);

    // ---- exhaustive: every (s, t) over a 2-letter alphabet, small lengths ----
    // two letters means substrings collide hard -- the walk sees long lcp runs
    // and lots of cross pairs, exactly where a half-check or off-by-one shows.
    {
        long long pairs = 0, diffs = 0;
        std::string first_s, first_t;
        for (int ls = 1; ls <= 7; ++ls) {
            const long long tots = 1LL << ls;  // 2^ls strings of length ls.
            for (long long cs = 0; cs < tots; ++cs) {
                std::string s(static_cast<std::size_t>(ls), 'a');
                for (int i = 0; i < ls; ++i)
                    s[static_cast<std::size_t>(i)] =
                        static_cast<char>('a' + ((cs >> i) & 1));
                for (int lt = 1; lt <= 7; ++lt) {
                    const long long tott = 1LL << lt;
                    for (long long ct = 0; ct < tott; ++ct) {
                        std::string t(static_cast<std::size_t>(lt), 'a');
                        for (int i = 0; i < lt; ++i)
                            t[static_cast<std::size_t>(i)] =
                                static_cast<char>('a' + ((ct >> i) & 1));
                        const int got = longest_common_substring(s, t);
                        const int want = reference_lcs(s, t);
                        ++pairs;
                        if (got != want) {
                            if (diffs == 0) {
                                first_s = s;
                                first_t = t;
                            }
                            ++diffs;
                        }
                    }
                }
            }
        }
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first diff: s=\"%s\" t=\"%s\"\n", first_s.c_str(),
                        first_t.c_str());
        }
        std::printf("exhaustive alpha2 |s|,|t| in 1..7: %lld pairs, %lld diffs\n",
                    pairs, diffs);
    }

    // ---- randomized differential vs the O(nm) dp, thousands of pairs ----
    // 2..4 letters, lengths up to ~60. small alphabets keep real common
    // substrings frequent, so the answer isn't trivially 0 or 1.
    {
        const int kCases = 40000;
        bool ok = cpf::differential(
            kCases, 0x1C5u,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 4);
                const std::int64_t ls = rng.in_range(1, 60);
                const std::int64_t lt = rng.in_range(1, 60);
                Pair p;
                p.s = make_string(ls, alpha, rng);
                p.t = make_string(lt, alpha, rng);
                return p;
            },
            [](const Pair& p) { return longest_common_substring(p.s, p.t); },
            [](const Pair& p) { return reference_lcs(p.s, p.t); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d cases, alpha 2..4, |s|,|t|<=60\n",
                    kCases);
    }

    // ---- planted-answer differential: force a long shared block ----
    // drop a random block into both strings so the answer is large, then let the
    // dp confirm the exact length. exercises the long-lcp path the random sweep
    // rarely reaches on its own.
    {
        const int kCases = 8000;
        bool ok = cpf::differential(
            kCases, 0xB10Cu,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 3);
                const std::int64_t blk = rng.in_range(1, 30);
                const std::string core = make_string(blk, alpha, rng);
                Pair p;
                p.s = make_string(rng.in_range(0, 15), alpha, rng) + core +
                      make_string(rng.in_range(0, 15), alpha, rng);
                p.t = make_string(rng.in_range(0, 15), alpha, rng) + core +
                      make_string(rng.in_range(0, 15), alpha, rng);
                return p;
            },
            [](const Pair& p) { return longest_common_substring(p.s, p.t); },
            [](const Pair& p) { return reference_lcs(p.s, p.t); });
        CPF_CHECK(ok);
        std::printf("planted-block differential: %d cases, shared core up to 30\n",
                    kCases);
    }

    return cpf::report();
}
