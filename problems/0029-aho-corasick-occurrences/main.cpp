#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

namespace {

// |T| up to 1e5 and sum|P| up to 1e5, all lowercase -- scanf per token would
// dominate the build. pull raw bytes through one buffer and hand back tokens.
struct FastIn {
    static constexpr int kBuf = 1 << 22;
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

    // one lowercase-letter token into out; false at EOF.
    bool word(std::string& out) {
        out.clear();
        int c = gc();
        while (c != -1 && (c < 'a' || c > 'z')) c = gc();
        if (c == -1) return false;
        while (c >= 'a' && c <= 'z') {
            out.push_back(static_cast<char>(c));
            c = gc();
        }
        return true;
    }

    bool integer(long long& out) {
        int c = gc();
        while (c != -1 && (c < '0' || c > '9')) c = gc();
        if (c == -1) return false;
        long long v = 0;
        while (c >= '0' && c <= '9') {
            v = v * 10 + (c - '0');
            c = gc();
        }
        out = v;
        return true;
    }
};

}  // namespace

// stdin: line 1 = T; line 2 = n; next n lines = patterns.
// stdout: the total occurrence count as one int64.
int main() {
    static FastIn in;
    std::string text;
    if (!in.word(text)) return 1;

    long long n = 0;
    if (!in.integer(n)) return 1;

    std::vector<std::string> patterns;
    patterns.reserve(static_cast<std::size_t>(n > 0 ? n : 0));
    std::string p;
    for (long long i = 0; i < n; ++i) {
        if (!in.word(p)) break;
        patterns.push_back(p);
    }

    const std::int64_t ans = p0029::total_occurrences(text, patterns);
    std::printf("%lld\n", static_cast<long long>(ans));
    return 0;
}
