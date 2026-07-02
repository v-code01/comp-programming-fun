#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// the third witness. not a dp -- the literal problem statement. pick any
// contiguous run of one letter, delete it, recurse, take the fewest operations.
// memoized on the string itself so it terminates on small inputs. it knows
// nothing about intervals, so agreeing with it proves the *recurrence*, not
// just that two dps copied each other.
std::unordered_map<std::string, int> g_brute_memo;

int brute(const std::string& s) {
    if (s.empty()) return 0;
    auto it = g_brute_memo.find(s);
    if (it != g_brute_memo.end()) return it->second;
    int best = static_cast<int>(s.size());  // deleting one letter at a time works.
    const int n = static_cast<int>(s.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            if (s[j] != s[i]) break;  // s[i..j] is a run of one letter.
            const int cand = 1 + brute(s.substr(0, i) + s.substr(j + 1));
            if (cand < best) best = cand;
        }
    }
    g_brute_memo[s] = best;
    return best;
}

// map a base-A integer onto a small lowercase alphabet -- cpf::Rng deals in
// ints, strings are what the problem eats.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    return s;
}

}  // namespace

int main() {
    using namespace p0006;
    ReferenceSolver ref;

    // the two statement examples -- ground truth.
    CPF_EQ(clear("abaca"), 3);
    CPF_EQ(clear("abcddcba"), 4);

    // edges, enumerated:
    CPF_EQ(clear("a"), 1);        // one letter.
    CPF_EQ(clear("aaaa"), 1);     // one run, one operation.
    CPF_EQ(clear("abcd"), 4);     // all distinct -- nothing merges, n operations.
    CPF_EQ(clear("ababab"), 4);   // two letters alternating, six long.
    CPF_EQ(clear("abcba"), 3);    // odd palindrome around a lone center.
    CPF_EQ(clear("abccba"), 3);   // even palindrome, matched middle.
    CPF_EQ(clear("abaa"), 2);     // first 'a' reaches the trailing run past 'b'.
    CPF_EQ(clear("aabaa"), 2);    // runs on both sides join across the center.

    // the counterexample that kills the tempting O(n^2). pairing only endpoints
    // (or folding into a neighbor) reports 3 here; the truth is 2, because the
    // leading 'a' has to reach a *middle* 'a'. keep it pinned.
    CPF_EQ(clear("abaa"), brute("abaa"));

    // exhaustive small: every string over a 3-letter alphabet up to length 8,
    // solution and reference both against the literal-deletion brute force. no
    // sampling -- 3^1..3^8 is 9840 strings and they all have to match.
    {
        long long total = 0, sol_diff = 0, ref_diff = 0;
        std::string s;
        for (int len = 1; len <= 8; ++len) {
            long long tot = 1;
            for (int i = 0; i < len; ++i) tot *= 3;
            for (long long code = 0; code < tot; ++code) {
                s.assign(static_cast<std::size_t>(len), 'a');
                long long c = code;
                for (int i = 0; i < len; ++i) {
                    s[static_cast<std::size_t>(i)] =
                        static_cast<char>('a' + (c % 3));
                    c /= 3;
                }
                const int b = brute(s);
                if (clear(s) != b) ++sol_diff;
                if (ref.solve(s) != b) ++ref_diff;
                ++total;
            }
        }
        CPF_CHECK(sol_diff == 0);
        CPF_CHECK(ref_diff == 0);
        std::printf("exhaustive alpha3 len1..8: %lld strings, sol_diff=%lld "
                    "ref_diff=%lld (vs deletion brute force)\n",
                    total, sol_diff, ref_diff);
    }

    // randomized differential: solution vs the last-char oracle, tens of
    // thousands of small strings over 2..4 letters, lengths 1..12. one gen
    // sweeps all three alphabets so no size is starved.
    const int kCases = 60000;
    bool ok = cpf::differential(
        kCases, 0xC1EA9,
        [](cpf::Rng& rng) {
            const std::int64_t alpha = rng.in_range(2, 4);
            const std::int64_t len = rng.in_range(1, 12);
            return make_string(len, alpha, rng);
        },
        [](const std::string& s) { return clear(s); },
        [&](const std::string& s) { return ref.solve(s); });
    CPF_CHECK(ok);
    std::printf("randomized differential: %d cases, alpha 2..4, len 1..12\n",
                kCases);

    return cpf::report();
}
