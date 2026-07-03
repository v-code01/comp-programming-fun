#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed integer scanner. scanf per token is too slow when n and sum k are
// both 1e5 -- that's millions of tokens. pull the whole stream once, parse ints
// by hand.
class Reader {
public:
    Reader() {
        len_ = std::fread(buf_, 1, sizeof(buf_), stdin);
    }

    // reads one non-negative integer. returns false at end of input.
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

// stdin: n, then n-1 edges, then q, then each query "k a_1 .. a_k".
// stdout: one line per query -- min captures, or -1.
int main() {
    Reader in;

    int n = 0;
    if (!in.next(n)) return 1;

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 0; i + 1 < n; ++i) {
        int u = 0, v = 0;
        in.next(u);
        in.next(v);
        edges.emplace_back(u, v);
    }

    p0016::Kingdom kingdom;
    kingdom.build(n, edges);

    int q = 0;
    in.next(q);

    // one output buffer, flushed once -- printf per line would dominate the wall.
    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 3 + 16);
    char tmp[16];

    std::vector<int> imp;
    for (int i = 0; i < q; ++i) {
        int k = 0;
        in.next(k);
        imp.clear();
        imp.resize(static_cast<std::size_t>(k));
        for (int j = 0; j < k; ++j) in.next(imp[static_cast<std::size_t>(j)]);

        const int ans = kingdom.query(imp);
        const int len = std::snprintf(tmp, sizeof(tmp), "%d\n", ans);
        out.append(tmp, static_cast<std::size_t>(len));
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
