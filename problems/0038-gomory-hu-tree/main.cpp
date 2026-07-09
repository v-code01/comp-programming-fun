#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <tuple>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner. up to n(n-1)/2 ~ 11k edges and 1e5 queries, so pull the
// stream in blocks and parse by hand -- scanf per token would dominate the wall.
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

// stdin:  line 1 "n m"; then m lines "u v c" (1-indexed, c >= 1); then a line
//         "q"; then q lines "u v".
// stdout: q lines -- the min cut between u and v.
int main() {
    Reader in;

    std::int64_t n64 = 0, m64 = 0;
    if (!in.next(n64) || !in.next(m64)) return 1;
    const int n = static_cast<int>(n64);
    const int m = static_cast<int>(m64);

    // sum parallel edges up front in a dense matrix -- n <= 150 makes it ~180 KB
    // and hands the solver one clean edge per unordered pair.
    std::vector<std::vector<std::int64_t>> cap(
        static_cast<std::size_t>(n),
        std::vector<std::int64_t>(static_cast<std::size_t>(n), 0));
    for (int e = 0; e < m; ++e) {
        std::int64_t u = 0, v = 0, c = 0;
        if (!in.next(u) || !in.next(v) || !in.next(c)) return 1;
        const int a = static_cast<int>(u) - 1;  // 1-indexed -> 0-indexed.
        const int b = static_cast<int>(v) - 1;
        if (a == b || a < 0 || b < 0 || a >= n || b >= n) continue;
        cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] += c;
        cap[static_cast<std::size_t>(b)][static_cast<std::size_t>(a)] += c;
    }

    std::vector<std::tuple<int, int, std::int64_t>> edges;
    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            if (cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] > 0) {
                edges.emplace_back(
                    a, b,
                    cap[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]);
            }
        }
    }

    const p0038::GomoryHu tree(n, edges);

    std::int64_t q64 = 0;
    if (!in.next(q64)) return 1;
    const int q = static_cast<int>(q64);

    // buffer the answers -- 1e5 printf calls flush far too often.
    std::vector<char> out;
    out.reserve(static_cast<std::size_t>(q) * 12);
    char tmp[24];
    for (int i = 0; i < q; ++i) {
        std::int64_t u = 0, v = 0;
        if (!in.next(u) || !in.next(v)) return 1;
        const std::int64_t ans =
            tree.min_cut(static_cast<int>(u) - 1, static_cast<int>(v) - 1);
        int len = std::snprintf(tmp, sizeof(tmp), "%lld\n",
                                static_cast<long long>(ans));
        out.insert(out.end(), tmp, tmp + len);
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
