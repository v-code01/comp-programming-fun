#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner. up to 5*n + 3m ~ 2e6 tokens, all non-negative, so pull
// the stream in blocks and parse digits by hand -- scanf per token would
// dominate the wall at this size.
class Reader {
public:
    Reader() { len_ = std::fread(buf_, 1, sizeof(buf_), stdin); }

    bool next(std::int64_t& out) {
        int c = get();
        while (c != EOF && (c < '0' || c > '9')) c = get();
        if (c == EOF) return false;
        std::int64_t v = 0;
        while (c >= '0' && c <= '9') {
            v = v * 10 + (c - '0');
            c = get();
        }
        out = v;
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

// stdin:  line 1 "n m s t k"; then m lines "u v w" for a directed edge u->v
//         with weight w >= 0.
// stdout: k lines, the k shortest s->t walk lengths non-decreasing, -1 for each
//         missing one.
int main() {
    Reader in;
    std::int64_t n, m, s, t, k;
    if (!in.next(n) || !in.next(m) || !in.next(s) || !in.next(t) || !in.next(k))
        return 1;

    std::vector<p0057::Edge> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (std::int64_t e = 0; e < m; ++e) {
        std::int64_t u, v, w;
        if (!in.next(u) || !in.next(v) || !in.next(w)) return 1;
        edges.push_back({static_cast<int>(u), static_cast<int>(v), w});
    }

    const std::vector<std::int64_t> ans = p0057::k_shortest_walks(
        static_cast<int>(n), static_cast<int>(s), static_cast<int>(t),
        static_cast<int>(k), edges);

    // one buffer, one write. printf per line would swamp the solve at k=1e5.
    std::string out;
    out.reserve(ans.size() * 8);
    char tmp[24];
    for (const std::int64_t x : ans) {
        const int len = std::snprintf(tmp, sizeof(tmp), "%lld",
                                      static_cast<long long>(x));
        out.append(tmp, static_cast<std::size_t>(len));
        out.push_back('\n');
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
