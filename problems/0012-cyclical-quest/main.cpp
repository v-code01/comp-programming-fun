#include <cstdint>
#include <cstdio>
#include <string>

#include "solution.hpp"

namespace {

// s is up to 1e6 chars, sum of query lengths up to 1e6, so scanf per token is
// too slow -- pull raw bytes through one big buffer and hand back whole tokens.
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

    // read one lowercase-letter token into out; false at EOF.
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

// stdin: line 1 = s, line 2 = q, next q lines = query strings.
// stdout: q lines, the position count cyclic-isomorphic to each query.
int main() {
    static FastIn in;
    std::string s;
    if (!in.word(s)) return 1;

    long long q = 0;
    if (!in.integer(q)) return 1;

    p0012::SuffixAutomaton sam;
    sam.build(s);

    std::string x;
    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 8);
    char tmp[24];
    for (long long i = 0; i < q; ++i) {
        if (!in.word(x)) break;
        const std::int64_t ans =
            sam.count_cyclic(x, static_cast<int>(i) + 1);
        const int m = std::snprintf(tmp, sizeof(tmp), "%lld\n",
                                    static_cast<long long>(ans));
        out.append(tmp, static_cast<std::size_t>(m));
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
