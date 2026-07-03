#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace {

// fread-backed scanner for signed integers. profits can be negative and there
// are up to n + 2m tokens (5000 + 40000), so pull the stream once and parse by
// hand -- scanf per token would dominate the wall.
class Reader {
public:
    Reader() { len_ = std::fread(buf_, 1, sizeof(buf_), stdin); }

    bool next(std::int64_t& out) {
        int c = get();
        while (c != EOF && c != '-' && (c < '0' || c > '9')) c = get();
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

// stdin:  line 1 "n m"; line 2 "p_1 .. p_n"; then m lines "i j" (select i -> need j).
// stdout: one int64 -- the max prerequisite-closed net profit.
int main() {
    Reader in;

    std::int64_t n64 = 0, m64 = 0;
    if (!in.next(n64) || !in.next(m64)) return 1;
    const int n = static_cast<int>(n64);
    const int m = static_cast<int>(m64);

    std::vector<std::int64_t> profits(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) in.next(profits[static_cast<std::size_t>(i)]);

    // ids arrive 1-indexed; the solver works 0-indexed.
    std::vector<std::pair<int, int>> prereqs;
    prereqs.reserve(static_cast<std::size_t>(m));
    for (int e = 0; e < m; ++e) {
        std::int64_t i = 0, j = 0;
        in.next(i);
        in.next(j);
        prereqs.emplace_back(static_cast<int>(i) - 1, static_cast<int>(j) - 1);
    }

    std::printf("%lld\n",
                static_cast<long long>(p0021::max_profit(profits, prereqs)));
    return 0;
}
