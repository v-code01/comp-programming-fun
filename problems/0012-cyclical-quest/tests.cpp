#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// cpf::Rng deals in ints; strings are what the automaton eats. map a base-alpha
// digit onto a lowercase alphabet -- small alphabets on purpose, so rotations
// collide and the dedupe path gets hammered.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    return s;
}

// one query against a freshly built automaton -- the differential's "solution".
std::int64_t sam_answer(const std::string& s, const std::string& x) {
    p0012::SuffixAutomaton sam;
    sam.build(s);
    return sam.count_cyclic(x, 1);
}

// enumerate every string over the first `alpha` letters with length in
// [1, max_len], calling f on each. base-alpha odometer, same trick as the
// exhaustive grid in 0006.
template <typename F>
void for_all_strings(int alpha, int max_len, F&& f) {
    std::string s;
    for (int len = 1; len <= max_len; ++len) {
        long long total = 1;
        for (int i = 0; i < len; ++i) total *= alpha;
        for (long long code = 0; code < total; ++code) {
            s.assign(static_cast<std::size_t>(len), 'a');
            long long c = code;
            for (int i = 0; i < len; ++i) {
                s[static_cast<std::size_t>(i)] =
                    static_cast<char>('a' + (c % alpha));
                c /= alpha;
            }
            f(s);
        }
    }
}

struct Case {
    std::string s;
    std::string x;
    bool operator==(const Case&) const = default;
};

}  // namespace

int main() {
    using namespace p0012;

    // ---- the statement example, verbatim. s = baabaabaaa, five queries. ----
    {
        const std::string s = "baabaabaaa";
        auto out = solve(s, {"a", "ba", "baa", "aabaa", "aaba"});
        CPF_EQ(out[0], std::int64_t{7});  // seven single 'a' positions.
        CPF_EQ(out[1], std::int64_t{5});  // ba/ab windows.
        CPF_EQ(out[2], std::int64_t{7});  // baa/aab/aba windows.
        CPF_EQ(out[3], std::int64_t{3});  // rotations of aabaa.
        CPF_EQ(out[4], std::int64_t{5});  // rotations of aaba.
        // the oracle must land on the same five -- ground truth, twice over.
        CPF_EQ(reference_count(s, "a"), std::int64_t{7});
        CPF_EQ(reference_count(s, "ba"), std::int64_t{5});
        CPF_EQ(reference_count(s, "baa"), std::int64_t{7});
        CPF_EQ(reference_count(s, "aabaa"), std::int64_t{3});
        CPF_EQ(reference_count(s, "aaba"), std::int64_t{5});
    }

    // ---- edges, enumerated. compare the SAM to the oracle so any drift shows. ----
    {
        const std::string s = "baabaabaaa";

        // single-char query -- just counts occurrences of that letter.
        CPF_EQ(sam_answer(s, "b"), reference_count(s, "b"));   // 3
        CPF_EQ(sam_answer(s, "c"), reference_count(s, "c"));   // 0, absent letter.

        // query longer than s -- no length-|x| window exists, answer 0.
        CPF_EQ(sam_answer("ab", "abc"), std::int64_t{0});
        CPF_EQ(sam_answer(s, "baabaabaaab"), std::int64_t{0});

        // all-same-letter s and x. every window matches; a run of five 'a's holds
        // six length-5 windows.
        CPF_EQ(sam_answer("aaaaa", "aaaaa"), std::int64_t{1});
        CPF_EQ(sam_answer("aaaaaaaaaa", "aaaaa"), std::int64_t{6});
        CPF_EQ(sam_answer("aaaaaaaaaa", "aaaaa"),
               reference_count("aaaaaaaaaa", "aaaaa"));

        // rotations that collide onto one state -- the dedupe has to fire.
        // "aa" has one distinct rotation; "aba" has three but two share letters.
        CPF_EQ(sam_answer("aaaa", "aa"), reference_count("aaaa", "aa"));   // 3
        CPF_EQ(sam_answer("ababab", "aba"),
               reference_count("ababab", "aba"));
        CPF_EQ(sam_answer("abababa", "aba"),
               reference_count("abababa", "aba"));

        // palindromic query -- rotations still distinct strings, still one state each.
        CPF_EQ(sam_answer("abacabaaba", "aba"),
               reference_count("abacabaaba", "aba"));
        CPF_EQ(sam_answer("racecarcar", "racecar"),
               reference_count("racecarcar", "racecar"));

        // query equal to s, and query longer only by rotation reach.
        CPF_EQ(sam_answer(s, s), reference_count(s, s));  // 1 (itself).
    }

    // ---- exhaustive: every (s, x) over a 2-letter alphabet, small lengths. ----
    // this is where dedupe/length-cap bugs surface: with only {a,b}, rotations
    // collapse constantly. no sampling -- every pair has to match the oracle.
    {
        long long pairs = 0, diffs = 0;
        std::string first_s, first_x;
        // s over lengths 1..8 (510 strings); for each build the SAM once and
        // sweep all queries x over lengths 1..6 (126 strings).
        for_all_strings(2, 8, [&](const std::string& s) {
            SuffixAutomaton sam;
            sam.build(s);
            int qid = 0;
            for_all_strings(2, 6, [&](const std::string& x) {
                ++qid;
                const std::int64_t got = sam.count_cyclic(x, qid);
                const std::int64_t want = reference_count(s, x);
                ++pairs;
                if (got != want && diffs == 0) {
                    first_s = s;
                    first_x = x;
                }
                if (got != want) ++diffs;
            });
        });
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first diff: s=\"%s\" x=\"%s\"\n", first_s.c_str(),
                        first_x.c_str());
        }
        std::printf("exhaustive alpha2  s.len1..8 x.len1..6: %lld pairs, %lld "
                    "diffs\n",
                    pairs, diffs);
    }

    // ---- randomized differential over 2..3 letters, longer strings. ----
    // 3-letter alphabet with lengths up to ~40 so states, clones, and repeated
    // rotations all get exercised at once. one gen sweeps both alphabets.
    {
        const int kCases = 50000;
        bool ok = cpf::differential(
            kCases, 0xCEEDu,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 3);
                const std::int64_t slen = rng.in_range(1, 40);
                const std::int64_t xlen = rng.in_range(1, 12);
                Case c;
                c.s = make_string(slen, alpha, rng);
                c.x = make_string(xlen, alpha, rng);
                return c;
            },
            [](const Case& c) { return sam_answer(c.s, c.x); },
            [](const Case& c) { return reference_count(c.s, c.x); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d cases, alpha 2..3, |s|<=40 "
                    "|x|<=12\n",
                    kCases);
    }

    // ---- multi-query dedupe across one automaton: ids must not bleed. ----
    // reuse one SAM for a batch of queries and check every answer still matches
    // the oracle -- proves visited_ tags are per-query, not sticky.
    {
        const std::string s = "abcabcabcabcxyzxyz";
        std::vector<std::string> qs = {"abc", "bca", "cab", "xyz", "yzx",
                                       "abcabc", "cabcab", "a", "z", "abcx"};
        auto got = solve(s, qs);
        long long diffs = 0;
        for (std::size_t i = 0; i < qs.size(); ++i) {
            if (got[i] != reference_count(s, qs[i])) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("multi-query one-automaton batch: %zu queries, %lld diffs\n",
                    qs.size(), diffs);
    }

    return cpf::report();
}
