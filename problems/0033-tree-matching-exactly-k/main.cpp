#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n K"; next n-1 lines = "u v w" (1-indexed tree edges, w >= 1).
// stdout: max weight of an exactly-K matching, or -1. answer is int64.
//
// n reaches 2e5 and ~6e5 integers arrive -- scanf per token is too slow, so pull
// the whole stream through one fread buffer and parse ints by hand. only w can
// exceed 32 bits on the tokens we build, but all fit in int64.
namespace {

struct Reader {
    static constexpr int kBuf = 1 << 20;
    char buf[kBuf];
    int len = 0, pos = 0;

    int getch() {
        if (pos == len) {
            len = static_cast<int>(std::fread(buf, 1, kBuf, stdin));
            pos = 0;
            if (len == 0) return -1;  // EOF.
        }
        return buf[pos++];
    }

    // next non-negative integer, or false at EOF. skips any separators.
    bool next(std::int64_t& out) {
        int c = getch();
        while (c != -1 && (c < '0' || c > '9')) c = getch();
        if (c == -1) return false;
        std::int64_t v = 0;
        while (c >= '0' && c <= '9') { v = v * 10 + (c - '0'); c = getch(); }
        out = v;
        return true;
    }
};

}  // namespace

int main() {
    Reader in;
    std::int64_t n64, k64;
    if (!in.next(n64) || !in.next(k64)) return 1;

    const int n = static_cast<int>(n64);
    std::vector<p0033::Edge> edges;
    if (n > 1) edges.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 0; i < n - 1; ++i) {
        std::int64_t u, v, w;
        if (!in.next(u) || !in.next(v) || !in.next(w)) return 1;
        edges.push_back({static_cast<int>(u), static_cast<int>(v), w});
    }

    std::printf("%lld\n",
                static_cast<long long>(p0033::solve(static_cast<int>(k64), n, edges)));
    return 0;
}
