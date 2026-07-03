#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace {

// n up to 1e6 and ~2e6 integers to read -- scanf per token is too slow, so pull
// raw bytes through one buffer and parse non-negative ints by hand.
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

    bool integer(int& out) {
        int c = gc();
        while (c != -1 && (c < '0' || c > '9')) c = gc();
        if (c == -1) return false;
        int v = 0;
        while (c >= '0' && c <= '9') {
            v = v * 10 + (c - '0');
            c = gc();
        }
        out = v;
        return true;
    }
};

}  // namespace

// stdin: line 1 = n; next n-1 lines = "u v" undirected edges, root = 1.
// stdout: n integers -- the dominant index of vertices 1..n, space separated.
int main() {
    static FastIn in;
    int n = 0;
    if (!in.integer(n)) return 1;

    std::vector<std::pair<int, int>> edges;
    if (n > 1) edges.reserve(static_cast<std::size_t>(n - 1));
    for (int i = 0; i < n - 1; ++i) {
        int u = 0, v = 0;
        if (!in.integer(u) || !in.integer(v)) return 1;
        edges.emplace_back(u, v);
    }

    const std::vector<int> ans = p0015::dominant_indices(n, edges);

    // one buffered write -- printf per answer would dominate the runtime.
    std::string out;
    out.reserve(ans.size() * 8);
    char tmp[16];
    for (std::size_t i = 0; i < ans.size(); ++i) {
        const int m = std::snprintf(tmp, sizeof(tmp), "%d%c", ans[i],
                                    i + 1 == ans.size() ? '\n' : ' ');
        out.append(tmp, static_cast<std::size_t>(m));
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
