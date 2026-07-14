#include <cstdio>
#include <string>

#include "solution.hpp"

// stdin: one line, the lowercase string s (up to 2e5 chars).
// stdout: distinct non-empty substrings, then the longest repeated length.
//
// scanf per character is wasted time at 2e5 -- one fread loop, keep the letters,
// drop the newline and anything else. answer (1) reaches ~2e10, so int64 the
// whole way out.
int main() {
    static char buf[1 << 16];
    std::string s;
    s.reserve(1 << 18);

    std::size_t got = 0;
    while ((got = std::fread(buf, 1, sizeof(buf), stdin)) > 0) {
        for (std::size_t i = 0; i < got; ++i) {
            const char c = buf[i];
            if (c >= 'a' && c <= 'z') s.push_back(c);
        }
    }

    const auto [distinct, longest] = p0040::suffix_tree_stats(s);
    std::printf("%lld %lld\n", static_cast<long long>(distinct),
                static_cast<long long>(longest));
    return 0;
}
