#include <cstdio>
#include <string>

#include "solution.hpp"

namespace {

// s and t run to 1e5 each, so read raw bytes through one buffer and hand back a
// lowercase token -- scanf per char would dominate the wall clock.
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

    // one run of lowercase letters into out; false at EOF.
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
};

}  // namespace

// stdin: line 1 = s, line 2 = t. stdout: length of their longest common
// contiguous substring, 0 if none.
int main() {
    static FastIn in;
    std::string s, t;
    if (!in.word(s)) return 1;
    if (!in.word(t)) return 1;
    std::printf("%d\n", p0028::longest_common_substring(s, t));
    return 0;
}
