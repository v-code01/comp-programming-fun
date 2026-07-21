#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: n, then n integers. stdout: the minimum L1 cost. n reaches 1e6, so a
// scanf per token is too slow -- pull the whole stream through one buffered
// reader and parse ints by hand.
namespace {
struct FastIn {
    static constexpr int kBuf = 1 << 16;
    char buf[kBuf];
    int len = 0, pos = 0;

    int gc() {
        if (pos == len) {
            len = static_cast<int>(std::fread(buf, 1, kBuf, stdin));
            pos = 0;
            if (len == 0) return -1;  // EOF.
        }
        return buf[pos++];
    }

    // read one signed int. returns false at clean EOF before any digit.
    bool read_int(int& out) {
        int c = gc();
        while (c != -1 && c != '-' && (c < '0' || c > '9')) c = gc();
        if (c == -1) return false;
        bool neg = false;
        if (c == '-') { neg = true; c = gc(); }
        int x = 0;
        while (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); c = gc(); }
        out = neg ? -x : x;
        return true;
    }
};
}  // namespace

int main() {
    FastIn in;
    int n;
    if (!in.read_int(n)) return 1;
    if (n < 0) return 1;

    std::vector<int> a(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!in.read_int(a[static_cast<std::size_t>(i)])) return 1;
    }

    std::printf("%lld\n", static_cast<long long>(p0050::solve(a)));
    return 0;
}
