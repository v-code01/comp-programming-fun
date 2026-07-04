#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: line 1 is "n K", line 2 is n prices. stdout: the max profit (int64).
// n reaches 1e6, so scanf per token is too slow -- pull the whole stream through
// one fread buffer and parse unsigned ints by hand. prices are >= 1, no sign.
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
    std::vector<int> prices(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        std::int64_t x;
        if (!in.next(x)) return 1;
        prices[static_cast<std::size_t>(i)] = static_cast<int>(x);
    }

    std::printf("%lld\n",
                static_cast<long long>(p0027::solve(static_cast<int>(k64), prices)));
    return 0;
}
