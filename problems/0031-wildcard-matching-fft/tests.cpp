#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// enumerate every string of length 1..max_len over the symbol set `syms`,
// calling f on each. base-|syms| odometer -- same trick 0012 uses. put '?' in
// syms for a pattern and every wildcard placement gets covered too.
template <typename F>
void for_all(const std::string& syms, int max_len, F&& f) {
    const int A = static_cast<int>(syms.size());
    std::string s;
    for (int len = 1; len <= max_len; ++len) {
        long long total = 1;
        for (int i = 0; i < len; ++i) total *= A;
        for (long long code = 0; code < total; ++code) {
            s.assign(static_cast<std::size_t>(len), syms[0]);
            long long c = code;
            for (int i = 0; i < len; ++i) {
                s[static_cast<std::size_t>(i)] =
                    syms[static_cast<std::size_t>(c % A)];
                c /= A;
            }
            f(s);
        }
    }
}

std::string rand_text(int n, int alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(n), 'a');
    for (auto& c : s)
        c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    return s;
}

std::string rand_pat(int m, int alpha, int wild_pct, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(m), 'a');
    for (auto& c : s) {
        if (rng.in_range(0, 99) < wild_pct)
            c = '?';
        else
            c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    }
    return s;
}

// the differential feeds both solvers one (text, pattern) pair.
struct Case {
    std::string text;
    std::string pat;
    bool operator==(const Case&) const = default;
};

}  // namespace

int main() {
    using namespace p0031;

    // ---- hand examples from the statement ----
    {
        // "a?c" over "abcabc": matches at 0 (abc) and 3 (abc).
        CPF_EQ(matches("abcabc", "a?c"), (std::vector<int>{0, 3}));

        // all-'?' pattern matches every window. |T|=6, |P|=3 -> starts 0..3.
        CPF_EQ(matches("abcabc", "???"), (std::vector<int>{0, 1, 2, 3}));

        // no match anywhere.
        CPF_EQ(matches("aaaaaa", "b"), (std::vector<int>{}));
        CPF_EQ(matches("abcabc", "xyz"), (std::vector<int>{}));

        // pattern equal to the whole text -> position 0 only.
        CPF_EQ(matches("abcabc", "abcabc"), (std::vector<int>{0}));

        // single-char pattern lists every occurrence of that letter.
        CPF_EQ(matches("abacaba", "a"), (std::vector<int>{0, 2, 4, 6}));

        // '?' at the ends, fixed middle.
        CPF_EQ(matches("xayxbyxcy", "?a?"), (std::vector<int>{0}));
        CPF_EQ(matches("banana", "?a?a?a"), (std::vector<int>{0}));

        // longer-than-text pattern: no window exists.
        CPF_EQ(matches("ab", "abc"), (std::vector<int>{}));

        // wildcard overlapping a run.
        CPF_EQ(matches("aaaa", "a?"), (std::vector<int>{0, 1, 2}));
    }

    // ---- oracle agreement on the hand cases, twice over ----
    {
        CPF_EQ(matches("abcabc", "a?c"), reference("abcabc", "a?c"));
        CPF_EQ(matches("mississippi", "?ss?"), reference("mississippi", "?ss?"));
        CPF_EQ(matches("mississippi", "ss"), reference("mississippi", "ss"));
    }

    // ---- exhaustive: every text over {a,b} and every pattern over {a,b,?} ----
    // small alphabet on purpose -- matches are common and wildcards land in
    // every slot. no sampling; each pair's full position vector must match.
    {
        long long pairs = 0, diffs = 0;
        std::string first_t, first_p;
        for_all("ab", 6, [&](const std::string& text) {
            for_all("ab?", 5, [&](const std::string& pat) {
                const std::vector<int> got = matches(text, pat);
                const std::vector<int> want = reference(text, pat);
                ++pairs;
                if (got != want) {
                    if (diffs == 0) {
                        first_t = text;
                        first_p = pat;
                    }
                    ++diffs;
                }
            });
        });
        CPF_CHECK(diffs == 0);
        if (diffs) {
            std::printf("  first diff: T=\"%s\" P=\"%s\"\n", first_t.c_str(),
                        first_p.c_str());
        }
        std::printf("exhaustive  T over {a,b} len1..6, P over {a,b,?} len1..5: "
                    "%lld pairs, %lld diffs\n",
                    pairs, diffs);
    }

    // ---- randomized differential, thousands of cases ----
    // alphabet 2..3 so windows collide, lengths up to ~80, ~25%% wildcards.
    // the harness replays any failing seed for you.
    {
        const int kCases = 30000;
        bool ok = cpf::differential(
            kCases, 0xC0FFEEu,  // arbitrary base seed.
            [](cpf::Rng& rng) {
                const int alpha = static_cast<int>(rng.in_range(2, 3));
                const int tlen = static_cast<int>(rng.in_range(1, 80));
                const int plen = static_cast<int>(rng.in_range(1, tlen));
                Case c;
                c.text = rand_text(tlen, alpha, rng);
                c.pat = rand_pat(plen, alpha, 25, rng);
                return c;
            },
            [](const Case& c) { return matches(c.text, c.pat); },
            [](const Case& c) { return reference(c.text, c.pat); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d cases, alpha 2..3, |T|<=80, "
                    "~25%% '?'\n",
                    kCases);
    }

    // ---- all-'?' patterns of every length against a random text ----
    // exercises the Csum=0 / weight-all-zero path: every window must match.
    {
        cpf::Rng rng(0xA11u);
        long long diffs = 0;
        for (int trial = 0; trial < 200; ++trial) {
            const int tlen = static_cast<int>(rng.in_range(1, 60));
            const int plen = static_cast<int>(rng.in_range(1, tlen));
            const std::string text = rand_text(tlen, 4, rng);
            const std::string pat(static_cast<std::size_t>(plen), '?');
            if (matches(text, pat) != reference(text, pat)) ++diffs;
        }
        CPF_CHECK(diffs == 0);
        std::printf("all-'?' patterns: 200 trials, %lld diffs\n", diffs);
    }

    // ---- large-scale precision: the FFT's actual danger zone ----
    // the differential above runs tiny transforms where rounding never bites.
    // these push |T|=2e5 through a full 2^18-2^19 transform at the WORST-case
    // magnitude (all 'z' -> every squared term is 26^2=676) and demand the exact
    // oracle answer. if the double FFT ever shaved a bit, it shows here.
    {
        const int N = 200000;

        // all-'z' text, all-'z' pattern -> every start matches. max magnitude.
        {
            const std::string text(static_cast<std::size_t>(N), 'z');
            const std::string pat(1000, 'z');
            CPF_CHECK(matches(text, pat) == reference(text, pat));
        }

        // same, but wildcards sprinkled in -- still every start matches.
        {
            const std::string text(static_cast<std::size_t>(N), 'z');
            std::string pat(1000, 'z');
            for (std::size_t i = 0; i < pat.size(); i += 7) pat[i] = '?';
            CPF_CHECK(matches(text, pat) == reference(text, pat));
        }

        // one 'a' breaks an otherwise all-'z' pattern -> zero matches.
        {
            const std::string text(static_cast<std::size_t>(N), 'z');
            std::string pat(1000, 'z');
            pat[500] = 'a';
            const std::vector<int> got = matches(text, pat);
            CPF_CHECK(got.empty());
            CPF_CHECK(got == reference(text, pat));
        }

        // random full-alphabet text, moderate pattern with wildcards, a few
        // seeds. oracle is O(n*m) ~ 4e8 per case -- fine for a handful.
        {
            long long diffs = 0;
            for (std::uint64_t seed = 1; seed <= 3; ++seed) {
                cpf::Rng rng(seed * 0x9E37u + 11u);
                const std::string text = rand_text(N, 26, rng);
                const std::string pat = rand_pat(2000, 26, 20, rng);
                if (matches(text, pat) != reference(text, pat)) ++diffs;
            }
            CPF_CHECK(diffs == 0);
            std::printf("large precision |T|=2e5: all-'z', wildcard, no-match, "
                        "3 random seeds -- %lld diffs\n",
                        diffs);
        }
    }

    return cpf::report();
}
