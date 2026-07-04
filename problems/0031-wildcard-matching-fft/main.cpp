#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

namespace {

// T and P are up to 2e5 chars each. one buffered read, hand back whole tokens --
// a token is a maximal run of lowercase letters or '?'. scanf per char would
// dominate the runtime the FFT is trying to save.
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

    static bool is_tok(int c) { return (c >= 'a' && c <= 'z') || c == '?'; }

    bool token(std::string& out) {
        out.clear();
        int c = gc();
        while (c != -1 && !is_tok(c)) c = gc();
        if (c == -1) return false;
        while (is_tok(c)) {
            out.push_back(static_cast<char>(c));
            c = gc();
        }
        return true;
    }
};

}  // namespace

// stdin:  line 1 = text T (lowercase), line 2 = pattern P (lowercase or '?').
// stdout: line 1 = count of matching 0-indexed start positions.
//         line 2 = those positions ascending, space-separated. empty line when
//                  the count is zero.
int main() {
    static FastIn in;
    std::string text, pat;
    if (!in.token(text)) return 1;
    if (!in.token(pat)) return 1;

    const std::vector<int> hits = p0031::matches(text, pat);

    std::string out;
    out.reserve(hits.size() * 7 + 16);
    char tmp[24];
    int m = std::snprintf(tmp, sizeof(tmp), "%zu\n", hits.size());
    out.append(tmp, static_cast<std::size_t>(m));
    for (std::size_t i = 0; i < hits.size(); ++i) {
        if (i) out.push_back(' ');
        m = std::snprintf(tmp, sizeof(tmp), "%d", hits[i]);
        out.append(tmp, static_cast<std::size_t>(m));
    }
    out.push_back('\n');  // always terminate line 2, even when empty.
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
