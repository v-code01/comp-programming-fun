#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

namespace {

// n, m up to 1e5, values up to ~1e9. buffered reader -- per-char scanf would
// swamp a run the tree keeps to a fraction of a second. unsigned decimals only,
// which is all the input holds.
struct FastIn {
    static constexpr int kBuf = 1 << 20;
    char buf[kBuf];
    int idx = 0;
    int len = 0;

    int gc() {
        if (idx >= len) {
            len = static_cast<int>(std::fread(buf, 1, kBuf, stdin));
            idx = 0;
            if (len <= 0) return -1;
        }
        return static_cast<unsigned char>(buf[idx++]);
    }

    bool read_i64(std::int64_t& out) {
        int c = gc();
        while (c != -1 && (c < '0' || c > '9')) c = gc();
        if (c == -1) return false;
        std::int64_t v = 0;
        while (c >= '0' && c <= '9') {
            v = v * 10 + (c - '0');
            c = gc();
        }
        out = v;
        return true;
    }
};

}  // namespace

// stdin:  line 1 = "n m". line 2 = a_0 … a_{n-1} in [0, p). line 3 = x_1 … x_m.
// stdout: A(x_1) … A(x_m) mod p, space-separated.
int main() {
    static FastIn in;
    std::int64_t n = 0, m = 0;
    if (!in.read_i64(n) || !in.read_i64(m) || n <= 0 || m <= 0) return 1;

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) {
        std::int64_t x = 0;
        if (!in.read_i64(x)) return 1;
        a[static_cast<std::size_t>(i)] = x % p0052::kMod;
    }

    std::vector<std::int64_t> pts(static_cast<std::size_t>(m));
    for (std::int64_t j = 0; j < m; ++j) {
        std::int64_t x = 0;
        if (!in.read_i64(x)) return 1;
        pts[static_cast<std::size_t>(j)] = x % p0052::kMod;
    }

    const std::vector<std::int64_t> res = p0052::multipoint_eval(a, pts);

    // one buffered write. up to 1e5 numbers under 1e9 -> ~10 chars each.
    std::string out;
    out.reserve(static_cast<std::size_t>(m) * 11 + 2);
    char tmp[24];
    for (std::int64_t j = 0; j < m; ++j) {
        if (j) out.push_back(' ');
        const int k = std::snprintf(
            tmp, sizeof(tmp), "%lld",
            static_cast<long long>(res[static_cast<std::size_t>(j)]));
        out.append(tmp, static_cast<std::size_t>(k));
    }
    out.push_back('\n');
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
