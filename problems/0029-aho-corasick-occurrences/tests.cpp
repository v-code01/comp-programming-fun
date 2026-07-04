#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// cpf::Rng deals in ints; the automaton eats strings. small alphabets on purpose
// -- with {a,b} patterns overlap and collide constantly, which is exactly where
// the fail links and the chain-count sum earn their keep or break.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s)
        c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    return s;
}

// one differential case: a text and its pattern list. output is a plain int64,
// so the harness diffs solution vs oracle directly.
struct Case {
    std::string text;
    std::vector<std::string> patterns;
};

}  // namespace

int main() {
    using namespace p0029;

    // ---- the statement example. abababab, patterns ab and aba. ----
    // "ab" starts at 0,2,4,6 -> 4. "aba" (len 3) at 0,2,4 -> 3. total 7.
    {
        const std::string t = "abababab";
        CPF_EQ(total_occurrences(t, {"ab", "aba"}), std::int64_t{7});
        // the oracle lands on the same 7 -- ground truth, twice over.
        CPF_EQ(reference_total(t, {"ab", "aba"}), std::int64_t{7});
        CPF_EQ(reference_total(t, {"ab"}), std::int64_t{4});
        CPF_EQ(reference_total(t, {"aba"}), std::int64_t{3});
    }

    // ---- absent pattern contributes nothing. ----
    CPF_EQ(total_occurrences("abababab", {"c"}), std::int64_t{0});
    CPF_EQ(total_occurrences("abababab", {"ba", "zzz"}),
           reference_total("abababab", {"ba", "zzz"}));  // 3 + 0.

    // ---- a duplicate pattern is counted once per listing. ----
    // "ab" x3 over abababab -> 4 each -> 12. the repeat must bump cnt_end again.
    CPF_EQ(total_occurrences("abababab", {"ab", "ab", "ab"}), std::int64_t{12});
    CPF_EQ(total_occurrences("aaaa", {"a", "a"}), std::int64_t{8});  // 4 + 4.

    // ---- a pattern equal to T occurs exactly once. ----
    CPF_EQ(total_occurrences("abcabc", {"abcabc"}), std::int64_t{1});
    CPF_EQ(total_occurrences("abcabc", {"abcabc", "abc"}),
           std::int64_t{1 + 2});  // whole + "abc" at 0 and 3.

    // ---- pattern longer than T -- no window, zero. ----
    CPF_EQ(total_occurrences("ab", {"abc"}), std::int64_t{0});

    // ---- one pattern is a suffix of another: the fail chain has to fire. ----
    // in "aXbXcabc" the state after "abc" fails through "bc"->"c", so a single
    // position must credit every pattern that ends up that chain at once.
    {
        const std::string t = "xabcx";
        auto pats = std::vector<std::string>{"abc", "bc", "c"};
        CPF_EQ(total_occurrences(t, pats), reference_total(t, pats));  // 1+1+1.
    }

    // ---- heavy overlap on one letter: 1e2 a's, several a-runs as patterns. ----
    {
        const std::string t(100, 'a');
        auto pats = std::vector<std::string>{"a", "aa", "aaa", "a"};
        CPF_EQ(total_occurrences(t, pats), reference_total(t, pats));
    }

    // ---- int64 guard: a large all-a text against a repeated single-a pattern.
    // 20000 a's x 300 copies of "a" = 6,000,000 -- past nothing dangerous yet,
    // but well past what an int-sized per-state sum would tolerate at scale. ----
    {
        const std::string t(20000, 'a');
        std::vector<std::string> pats(300, "a");
        CPF_EQ(total_occurrences(t, pats), std::int64_t{20000} * 300);
    }

    // ---- randomized differential: thousands of small cases. ----
    // alphabet 2..3, short text, several short patterns with repeats allowed so
    // overlaps, duplicates and absent patterns all show up. every case must match
    // the naive oracle; the first disagreement prints its seed to replay.
    {
        const int kCases = 40000;
        bool ok = cpf::differential(
            kCases, 0xA0C0u,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 3);
                const std::int64_t tlen = rng.in_range(1, 25);
                const std::int64_t n = rng.in_range(1, 8);
                Case c;
                c.text = make_string(tlen, alpha, rng);
                c.patterns.reserve(static_cast<std::size_t>(n));
                for (std::int64_t i = 0; i < n; ++i) {
                    const std::int64_t plen = rng.in_range(1, 6);
                    c.patterns.push_back(make_string(plen, alpha, rng));
                }
                return c;
            },
            [](const Case& c) { return total_occurrences(c.text, c.patterns); },
            [](const Case& c) { return reference_total(c.text, c.patterns); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d cases, alpha 2..3, |T|<=25, "
                    "n<=8, |P|<=6\n",
                    kCases);
    }

    // ---- second differential, single-letter alphabet: maximal overlap. ----
    // everything is 'a', so every pattern is a run and every state is on one
    // fail chain -- the count-sum has the most work and the most chances to
    // double-count. this is where an off-by-one in sum_end would light up.
    {
        const int kCases = 10000;
        bool ok = cpf::differential(
            kCases, 0xBEE5u,
            [](cpf::Rng& rng) {
                const std::int64_t tlen = rng.in_range(1, 40);
                const std::int64_t n = rng.in_range(1, 10);
                Case c;
                c.text = make_string(tlen, 1, rng);  // all 'a'.
                c.patterns.reserve(static_cast<std::size_t>(n));
                for (std::int64_t i = 0; i < n; ++i)
                    c.patterns.push_back(
                        make_string(rng.in_range(1, 8), 1, rng));
                return c;
            },
            [](const Case& c) { return total_occurrences(c.text, c.patterns); },
            [](const Case& c) { return reference_total(c.text, c.patterns); });
        CPF_CHECK(ok);
        std::printf("single-letter differential: %d cases, |T|<=40, n<=10, "
                    "|P|<=8\n",
                    kCases);
    }

    return cpf::report();
}
