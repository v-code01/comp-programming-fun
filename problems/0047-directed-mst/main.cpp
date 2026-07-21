#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner. up to n + 3m ~ 3e5 tokens, weights up to 1e9, so pull
// the stream in blocks and parse digits by hand -- scanf per token would
// dominate the wall at this size. reads into int64 so a weight never truncates.
class Reader {
public:
    Reader() { len_ = std::fread(buf_, 1, sizeof(buf_), stdin); }

    bool next(std::int64_t& out) {
        int c = get();
        while (c != EOF && (c < '0' || c > '9') && c != '-') c = get();
        if (c == EOF) return false;
        bool neg = false;
        if (c == '-') {
            neg = true;
            c = get();
        }
        std::int64_t v = 0;
        while (c >= '0' && c <= '9') {
            v = v * 10 + (c - '0');
            c = get();
        }
        out = neg ? -v : v;
        return true;
    }

private:
    int get() {
        if (pos_ >= len_) {
            len_ = std::fread(buf_, 1, sizeof(buf_), stdin);
            pos_ = 0;
            if (len_ == 0) return EOF;
        }
        return buf_[pos_++];
    }

    char buf_[1 << 16];
    std::size_t len_ = 0;
    std::size_t pos_ = 0;
};

}  // namespace

// stdin:  line 1 "n m r"; then m lines "u v w" for a directed edge u->v.
// stdout: the minimum arborescence weight rooted at r, or -1 if none exists.
int main() {
    Reader in;
    std::int64_t n = 0, m = 0, r = 0;
    if (!in.next(n) || !in.next(m) || !in.next(r)) return 1;

    std::vector<p0047::Edge> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (std::int64_t i = 0; i < m; ++i) {
        std::int64_t u = 0, v = 0, w = 0;
        if (!in.next(u) || !in.next(v) || !in.next(w)) return 1;
        edges.push_back(
            p0047::Edge{static_cast<int>(u), static_cast<int>(v), w});
    }

    const std::int64_t ans =
        p0047::min_arborescence(static_cast<int>(n), edges, static_cast<int>(r));
    std::printf("%lld\n", static_cast<long long>(ans));
    return 0;
}
