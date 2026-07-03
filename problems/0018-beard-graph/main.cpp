#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = n; next n-1 lines = edges "u v" (1-based, edge i is line i);
// next line = m; then m queries -- "1 i" paints edge i black, "2 i" paints it
// white, "3 a b" prints the black-path edge count from a to b (or -1, or 0 when
// a == b). stdout: one line per type-3 query.
//
// n, m up to 1e5. scanf per token drags at that many tokens, so pull raw bytes
// through one buffer and stage all output in one string flushed once. edges start
// black -- build() folds that in.
namespace {

constexpr int kBuf = 1 << 16;
char g_in[kBuf];
int g_pos = 0;
int g_len = 0;

inline int gc() {
    if (g_pos == g_len) {
        g_len = static_cast<int>(std::fread(g_in, 1, kBuf, stdin));
        g_pos = 0;
        if (g_len == 0) return -1;
    }
    return g_in[g_pos++];
}

// every integer here is a positive vertex/edge index or a type in {1,2,3}.
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
    const int n = rd();

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 0; i < n - 1; ++i) {
        const int u = rd();
        const int v = rd();
        edges.emplace_back(u - 1, v - 1);  // to 0-based.
    }

    p0018::BeardHLD tree;
    tree.build(n, edges);

    const int m = rd();
    std::string out;
    out.reserve(static_cast<std::size_t>(m) * 3);  // rough, dodges regrowth.
    char tmp[16];

    for (int i = 0; i < m; ++i) {
        const int t = rd();
        if (t == 1) {
            tree.paint_black(rd() - 1);  // edge index to 0-based.
        } else if (t == 2) {
            tree.paint_white(rd() - 1);
        } else {
            const int a = rd() - 1;  // vertices to 0-based.
            const int b = rd() - 1;
            const int len =
                std::snprintf(tmp, sizeof(tmp), "%d\n", tree.query(a, b));
            out.append(tmp, static_cast<std::size_t>(len));
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
