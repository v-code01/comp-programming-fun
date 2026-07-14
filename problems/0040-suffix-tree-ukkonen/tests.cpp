#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using Answer = std::pair<std::int64_t, std::int64_t>;

Answer got(const std::string& s) { return p0040::suffix_tree_stats(s); }
Answer want(const std::string& s) { return p0040::reference_stats(s); }

Answer pair_of(std::int64_t distinct, std::int64_t longest) {
    return {distinct, longest};
}

// cpf::Rng deals in ints; the tree eats letters. small alphabets on purpose --
// with 2 or 3 letters the repeats are dense, edges split constantly, and the
// suffix-link hops actually get exercised. 26 random letters is the easy case.
std::string make_string(std::int64_t len, std::int64_t alpha, cpf::Rng& rng) {
    std::string s(static_cast<std::size_t>(len), 'a');
    for (auto& c : s) {
        c = static_cast<char>('a' + rng.in_range(0, alpha - 1));
    }
    return s;
}

// every string over the first `alpha` letters with length in [1, max_len].
// a base-alpha odometer.
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
    // ---- hand examples, checked on paper first. ----
    // banana: len1 {b,a,n}=3, len2 {ba,an,na}=3, len3 {ban,ana,nan}=3,
    // len4 {bana,anan,nana}=3, len5 {banan,anana}=2, len6 {banana}=1 -> 15.
    // "ana" sits at 1 and at 3 -- overlapping, still two occurrences -> 3.
    CPF_EQ(got("banana"), pair_of(15, 3));
    CPF_EQ(want("banana"), pair_of(15, 3));

    // aaaa: distinct is one per length -- a, aa, aaa, aaaa. "aaa" occurs at 0
    // and at 1 -> 3. the all-overlap case.
    CPF_EQ(got("aaaa"), pair_of(4, 3));
    CPF_EQ(want("aaaa"), pair_of(4, 3));

    // abc: every substring unique -- 6 of them, nothing repeats.
    CPF_EQ(got("abc"), pair_of(6, 0));
    CPF_EQ(want("abc"), pair_of(6, 0));

    // ---- edges the tree has to survive. ----
    CPF_EQ(got("a"), pair_of(1, 0));         // single character. no internal node.
    CPF_EQ(got("aa"), pair_of(2, 1));        // shortest possible repeat.
    CPF_EQ(got("ab"), pair_of(3, 0));        // two distinct letters, no repeat.
    CPF_EQ(got("abababab"), pair_of(15, 6)); // periodic: "ababab" repeats.
    CPF_EQ(got("abcdefghijklmnopqrstuvwxyz"),
           pair_of(26 * 27 / 2, 0));         // all 26 letters, all distinct.
    CPF_EQ(got("aabaacaad"), want("aabaacaad"));   // shared prefixes, splits.
    CPF_EQ(got("abcabcabc"), want("abcabcabc"));   // period 3.
    CPF_EQ(got("aaaaaaaaaa"), pair_of(10, 9));     // mono, longer.
    CPF_EQ(got(std::string(64, 'z')), pair_of(64, 63));

    // the oracle has to agree on the shapes too -- ground truth, twice over.
    CPF_EQ(want("abababab"), pair_of(15, 6));
    CPF_EQ(want("abcdefghijklmnopqrstuvwxyz"), pair_of(26 * 27 / 2, 0));

    // ---- exhaustive: every 2-letter string, lengths 1..14. ----
    // {a,b} only, so suffixes collide constantly -- exactly where a botched
    // active-point reset or a missing suffix link shows up.
    {
        long long strings = 0, diffs = 0;
        std::string first;
        for_all_strings(2, 14, [&](const std::string& s) {
            ++strings;
            if (got(s) != want(s)) {
                if (diffs == 0) first = s;
                ++diffs;
            }
        });
        CPF_CHECK(diffs == 0);
        if (diffs) std::printf("  first diff: s=\"%s\"\n", first.c_str());
        std::printf("exhaustive alpha2 len1..14: %lld strings, %lld diffs\n",
                    strings, diffs);
    }

    // ---- exhaustive: every 3-letter string, lengths 1..9. ----
    {
        long long strings = 0, diffs = 0;
        std::string first;
        for_all_strings(3, 9, [&](const std::string& s) {
            ++strings;
            if (got(s) != want(s)) {
                if (diffs == 0) first = s;
                ++diffs;
            }
        });
        CPF_CHECK(diffs == 0);
        if (diffs) std::printf("  first diff: s=\"%s\"\n", first.c_str());
        std::printf("exhaustive alpha3 len1..9: %lld strings, %lld diffs\n",
                    strings, diffs);
    }

    // ---- randomized differential, longer strings over 2..4 letters. ----
    // up to 120 characters: long enough that the tree is deep, the active point
    // walks down through several edges, and the link hops telescope many times
    // over. that amortized walk is where ukkonen dies if it's wrong.
    {
        const int kCases = 20000;
        bool ok = cpf::differential(
            kCases, 0x40FEu,
            [](cpf::Rng& rng) {
                const std::int64_t alpha = rng.in_range(2, 4);
                const std::int64_t len = rng.in_range(1, 120);
                return make_string(len, alpha, rng);
            },
            [](const std::string& s) { return got(s); },
            [](const std::string& s) { return want(s); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d strings, alpha 2..4, |s|<=120\n",
                    kCases);
    }

    // ---- randomized differential, the wide alphabet. ----
    // 26 letters: repeats are rare, most extensions are plain rule-2 leaves, and
    // the deepest internal node is shallow. the other half of the state space.
    {
        const int kCases = 5000;
        bool ok = cpf::differential(
            kCases, 0x1A26u,
            [](cpf::Rng& rng) {
                const std::int64_t len = rng.in_range(1, 120);
                return make_string(len, 26, rng);
            },
            [](const std::string& s) { return got(s); },
            [](const std::string& s) { return want(s); });
        CPF_CHECK(ok);
        std::printf("randomized differential: %d strings, alpha 26, |s|<=120\n",
                    kCases);
    }

    // ---- structured shapes, where random strings never go. ----
    // all-same, periodic at every period, a fibonacci word (worst case for
    // repeat structure), and a thue-morse-ish prefix. each against the oracle.
    {
        long long shapes = 0, diffs = 0;
        auto check = [&](const std::string& s) {
            ++shapes;
            if (got(s) != want(s)) ++diffs;
        };

        for (int n = 1; n <= 60; ++n) check(std::string(static_cast<std::size_t>(n), 'a'));

        for (int period = 1; period <= 8; ++period) {
            for (int n = 1; n <= 60; ++n) {
                std::string s;
                for (int i = 0; i < n; ++i) {
                    s.push_back(static_cast<char>('a' + i % period));
                }
                check(s);
            }
        }

        // fibonacci word: f(k) = f(k-1) + f(k-2). packed with repeats at every
        // scale, and no period at all.
        {
            std::string a = "a", b = "ab";
            for (int k = 0; k < 8; ++k) {
                check(b);
                std::string next = b + a;
                a = b;
                b = next;
            }
        }

        // thue-morse: no cube, so the deepest internal node stays shallow while
        // the substring count runs high. the opposite stress from all-'a'.
        {
            std::string tm = "a";
            for (int k = 0; k < 8; ++k) {
                std::string flip = tm;
                for (auto& c : flip) c = (c == 'a') ? 'b' : 'a';
                tm += flip;
                if (tm.size() <= 128) check(tm);
            }
        }

        // every suffix of the alphabet -- all-distinct, so answer (2) is 0 and
        // answer (1) is n(n+1)/2 exactly.
        for (int n = 1; n <= 26; ++n) {
            check(std::string("abcdefghijklmnopqrstuvwxyz").substr(
                0, static_cast<std::size_t>(n)));
        }

        CPF_CHECK(diffs == 0);
        std::printf("structured shapes: %lld strings, %lld diffs\n", shapes,
                    diffs);
    }

    return cpf::report();
}
