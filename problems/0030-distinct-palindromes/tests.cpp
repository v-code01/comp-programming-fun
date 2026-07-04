#include <cstdint>
#include <cstdio>
#include <string>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

// cpf::Rng hands back ints; the eertree eats strings. map each int onto the first
// `alpha` lowercase letters. keep alpha small -- 2 or 3 -- so palindromes get
// dense and the suffix-link walk and the dedupe both get hammered.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s) {
        c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    }
    return s;
}

// enumerate every string over the first `alpha` letters with length in
// [1, max_len], calling f on each. a base-alpha odometer -- same trick as 0012.
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

}  // namespace

int main() {
    using namespace p0030;

    // ---- hand examples from the statement. ----
    CPF_EQ(count_distinct("aba"), std::int64_t{3});   // a, b, aba.
    CPF_EQ(count_distinct("aabb"), std::int64_t{4});  // a, aa, b, bb.
    CPF_EQ(count_distinct("a"), std::int64_t{1});     // a.
    CPF_EQ(count_distinct("aaaa"), std::int64_t{4});  // a, aa, aaa, aaaa.

    // the oracle has to land on the same four -- ground truth, twice over.
    CPF_EQ(reference_count("aba"), std::int64_t{3});
    CPF_EQ(reference_count("aabb"), std::int64_t{4});
    CPF_EQ(reference_count("a"), std::int64_t{1});
    CPF_EQ(reference_count("aaaa"), std::int64_t{4});

    // ---- edges the tree has to get right. ----
    CPF_EQ(count_distinct(""), std::int64_t{0});   // empty -- just the two roots.
    CPF_EQ(count_distinct("ab"), std::int64_t{2}); // a, b -- no even palindrome.
    CPF_EQ(count_distinct("abc"), std::int64_t{3});          // three singletons.
    CPF_EQ(count_distinct("abba"), std::int64_t{4});         // a, b, bb, abba.
    CPF_EQ(count_distinct("racecar"), std::int64_t{7});      // full odd nest.
    CPF_EQ(count_distinct("aaaaa"), std::int64_t{5});        // one per length.

    // each of those against the oracle too, so any drift shows.
    CPF_EQ(count_distinct("abba"), reference_count("abba"));
    CPF_EQ(count_distinct("racecar"), reference_count("racecar"));
    CPF_EQ(count_distinct("aaaaa"), reference_count("aaaaa"));

    // ---- exhaustive: every string over a 2-letter alphabet, lengths 1..16. ----
    // with only {a,b} palindromes collide constantly -- exactly where a bad
    // suffix-link walk or a mishandled odd root would show. no sampling.
    {
        long long strings = 0, diffs = 0;
        std::string first;
        for_all_strings(2, 16, [&](const std::string& s) {
            ++strings;
            if (count_distinct(s) != reference_count(s)) {
                if (diffs == 0) first = s;
                ++diffs;
            }
        });
        CPF_CHECK(diffs == 0);
        if (diffs) std::printf("  first diff: s=\"%s\"\n", first.c_str());
        std::printf("exhaustive alpha2 len1..16: %lld strings, %lld diffs\n",
                    strings, diffs);
    }

    // ---- exhaustive: every string over a 3-letter alphabet, lengths 1..10. ----
    {
        long long strings = 0, diffs = 0;
        std::string first;
        for_all_strings(3, 10, [&](const std::string& s) {
            ++strings;
            if (count_distinct(s) != reference_count(s)) {
                if (diffs == 0) first = s;
                ++diffs;
            }
        });
        CPF_CHECK(diffs == 0);
        if (diffs) std::printf("  first diff: s=\"%s\"\n", first.c_str());
        std::printf("exhaustive alpha3 len1..10: %lld strings, %lld diffs\n",
                    strings, diffs);
    }

    // ---- randomized differential, longer strings over 2..3 letters. ----
    // lengths up to 80 so the link walk climbs deep and telescopes back down many
    // times over -- the amortization is where the bug-prone walk earns its keep.
    {
        const int kCases = 40000;
        bool ok = cpf::differential(
            kCases, 0xEE27u,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 3);
                const std::int64_t slen = rng.in_range(1, 80);
                return make_string(slen, alpha, rng);
            },
            [](const std::string& s) { return count_distinct(s); },
            [](const std::string& s) { return reference_count(s); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d strings, alpha 2..3, |s|<=80\n",
                    kCases);
    }

    return cpf::report();
}
