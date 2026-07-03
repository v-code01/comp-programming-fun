#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

namespace {

// n colors plus 2*(n-1) edge endpoints is ~3e5 integers at the cap -- pull them
// through one buffered reader so parsing never shows up next to the O(n log n)
// solve.
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

// stdin: line 1 = n; line 2 = n colors; next n-1 lines = undirected edges "u v",
// root is vertex 1. stdout: n answers for vertices 1..n, space-separated.
int main() {
    static FastIn in;
    int n = 0;
    if (!in.integer(n)) return 1;

    std::vector<int> color(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        if (!in.integer(color[static_cast<std::size_t>(i)])) return 1;
    }

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 0; i < n - 1; ++i) {
        int u = 0, v = 0;
        if (!in.integer(u) || !in.integer(v)) return 1;
        edges.emplace_back(u, v);
    }

    std::vector<std::int64_t> ans = p0013::solve(n, color, edges);

    // one string, one fwrite -- n separate printf calls would dominate the runtime.
    std::string out;
    out.reserve(static_cast<std::size_t>(n) * 8);
    char tmp[24];
    for (int i = 0; i < n; ++i) {
        int m = std::snprintf(tmp, sizeof(tmp), "%lld%c",
                              static_cast<long long>(ans[static_cast<std::size_t>(i)]),
                              i + 1 < n ? ' ' : '\n');
        out.append(tmp, static_cast<std::size_t>(m));
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
