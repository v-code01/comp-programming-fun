#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner. up to n + 2m ~ 2e5 + 6e5 tokens, all non-negative, so
// pull the stream in blocks and parse digits by hand -- scanf per token would
// dominate the wall at this size.
class Reader {
public:
    Reader() { len_ = std::fread(buf_, 1, sizeof(buf_), stdin); }

    bool next(int& out) {
        int c = get();
        while (c != EOF && (c < '0' || c > '9')) c = get();
        if (c == EOF) return false;
        int v = 0;
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

// stdin:  line 1 "n m"; then m lines "u v" for a directed edge u->v.
// stdout: idom(1..n), space-separated. 0 for the root and unreachable vertices.
int main() {
    Reader in;
    int n = 0, m = 0;
    if (!in.next(n) || !in.next(m)) return 1;

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(m));
    for (int e = 0; e < m; ++e) {
        int u = 0, v = 0;
        if (!in.next(u) || !in.next(v)) return 1;
        edges.emplace_back(u, v);
    }

    const std::vector<int> idom = p0036::dominator_tree(n, edges, 1);

    // one big buffer, one write. printf per integer would swamp the solve.
    std::string out;
    out.reserve(static_cast<std::size_t>(n) * 7);
    char tmp[12];
    for (int v = 1; v <= n; ++v) {
        int len = std::snprintf(tmp, sizeof(tmp), "%d", idom[static_cast<std::size_t>(v)]);
        out.append(tmp, static_cast<std::size_t>(len));
        out.push_back(v == n ? '\n' : ' ');
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
