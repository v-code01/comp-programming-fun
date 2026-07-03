#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

// stdin: "n m", then n colors (1-indexed vertices 1..n), then n-1 edges "u v",
// then m queries --
//   1 v c  -> repaint every vertex in v's subtree the color c
//   2 v    -> print the number of distinct colors in v's subtree
// one line of output per type-2 query.
//
// fast IO. n, m up to 4e5 means ~1e6 integers to read; scanf per token would
// dominate the run, so slurp stdin through a byte buffer and stage all output in
// one string flushed once at the end.
namespace {

constexpr int kBuf = 1 << 16;
char g_in[kBuf];
int g_pos = 0, g_len = 0;

inline int gc() {
    if (g_pos == g_len) {
        g_len = static_cast<int>(std::fread(g_in, 1, kBuf, stdin));
        g_pos = 0;
        if (g_len == 0) return -1;  // EOF.
    }
    return g_in[g_pos++];
}

// every integer in this problem is non-negative.
inline int rd() {
    int c = gc();
    while (c != -1 && (c < '0' || c > '9')) c = gc();
    int x = 0;
    while (c >= '0' && c <= '9') {
        x = x * 10 + (c - '0');
        c = gc();
    }
    return x;
}

}  // namespace

int main() {
    int n = rd();
    int m = rd();

    std::vector<int> color(static_cast<std::size_t>(n) + 1, 0);
    for (int i = 1; i <= n; ++i) color[static_cast<std::size_t>(i)] = rd();

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 0; i < n - 1; ++i) {
        int u = rd();
        int v = rd();
        edges.emplace_back(u, v);
    }

    p0017::Tree tree;
    tree.build(n, color, edges);

    std::string out;
    out.reserve(static_cast<std::size_t>(m) * 3);  // rough, avoids regrowth.

    for (int q = 0; q < m; ++q) {
        int type = rd();
        if (type == 1) {
            int v = rd();
            int c = rd();
            tree.repaint(v, c);
        } else {
            int v = rd();
            out += std::to_string(tree.distinct(v));
            out += '\n';
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
