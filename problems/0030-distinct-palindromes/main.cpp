#include <cstdio>
#include <string>

#include "solution.hpp"

// stdin: one line, the lowercase string s (up to 1e6 chars). stdout: the count of
// distinct palindromic substrings.
//
// scanf per char is too slow at 1e6 -- pull the whole input through one fread
// buffer and keep only the letters.
int main() {
    static char buf[1 << 22];
    std::string s;
    s.reserve(1 << 20);

    std::size_t got = 0;
    while ((got = std::fread(buf, 1, sizeof(buf), stdin)) > 0) {
        for (std::size_t i = 0; i < got; ++i) {
            const char c = buf[i];
            if (c >= 'a' && c <= 'z') s.push_back(c);
        }
    }

    std::printf("%lld\n", static_cast<long long>(p0030::count_distinct(s)));
    return 0;
}
