#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner. up to n(n-1)/2 ~ 125k edges, three tokens each, so pull
// the stream in blocks and parse by hand -- scanf per token would dominate.
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

// stdin:  line 1 "n m"; then m lines "u v w" (1-indexed vertices, w >= 1).
// stdout: one int64 -- the maximum total weight of a matching (need not cover
// every vertex; empty matching = 0).
int main() {
    Reader in;

    std::int64_t n64 = 0, m64 = 0;
    if (!in.next(n64) || !in.next(m64)) return 1;
    const int n = static_cast<int>(n64);
    const int m = static_cast<int>(m64);

    // dedup parallel edges by max up front, in a dense matrix -- n <= 500 makes
    // this 2 MB and keeps the solver's input clean of duplicates.
    std::vector<std::vector<std::int64_t>> w(
        static_cast<std::size_t>(n),
        std::vector<std::int64_t>(static_cast<std::size_t>(n), 0));
    for (int e = 0; e < m; ++e) {
        std::int64_t u = 0, v = 0, ww = 0;
        if (!in.next(u) || !in.next(v) || !in.next(ww)) return 1;
        const int a = static_cast<int>(u) - 1;  // 1-indexed -> 0-indexed.
        const int b = static_cast<int>(v) - 1;
        if (a == b || a < 0 || b < 0 || a >= n || b >= n) continue;
        if (ww > w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]) {
            w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = ww;
            w[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] = ww;
        }
    }

    std::vector<std::tuple<int, int, std::int64_t>> edges;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            if (w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] > 0)
                edges.emplace_back(
                    a, b, w[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]);

    std::printf("%lld\n",
                static_cast<long long>(p0032::max_weight_matching(n, edges)));
    return 0;
}
