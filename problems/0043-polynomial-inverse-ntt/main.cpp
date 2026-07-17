#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

namespace {

// n up to 1e5, coefficients up to ~1e9. buffered reader -- per-char scanf would
// dominate a run the NTT keeps to milliseconds. reads unsigned decimals only,
// which is all the input holds (n, then non-negative coefficients).
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

// stdin:  line 1 = n. line 2 = a_0 … a_{n-1}, each in [0, p), a_0 != 0.
// stdout: b_0 … b_{n-1}, space-separated, where B = A^{-1} mod x^n.
int main() {
    static FastIn in;
    std::int64_t n = 0;
    if (!in.read_i64(n) || n <= 0) return 1;

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) {
        std::int64_t x = 0;
        if (!in.read_i64(x)) return 1;
        a[static_cast<std::size_t>(i)] = x % p0043::kMod;
    }

    const std::vector<std::int64_t> b =
        p0043::inverse(a, static_cast<int>(n));

    // one buffered write. up to 1e5 numbers under 1e9 -> ~10 chars each.
    std::string out;
    out.reserve(static_cast<std::size_t>(n) * 11 + 2);
    char tmp[24];
    for (std::int64_t i = 0; i < n; ++i) {
        if (i) out.push_back(' ');
        const int m = std::snprintf(
            tmp, sizeof(tmp), "%lld",
            static_cast<long long>(b[static_cast<std::size_t>(i)]));
        out.append(tmp, static_cast<std::size_t>(m));
    }
    out.push_back('\n');
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
