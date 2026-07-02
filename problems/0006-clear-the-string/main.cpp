#include <cstdio>
#include <string>

#include "solution.hpp"

// stdin: line 1 is n, line 2 is the string. stdout: min operations to clear it.
// n is read and discarded -- the string carries its own length, and trusting
// the count over the bytes is how off-by-ones sneak in.
int main() {
    int n;
    if (std::scanf("%d", &n) != 1) return 1;
    std::string s;
    // %s stops at whitespace, which is exactly the token we want on line 2.
    char buf[512];
    if (std::scanf("%511s", buf) != 1) return 1;
    s = buf;
    std::printf("%d\n", p0006::clear(s));
    return 0;
}
