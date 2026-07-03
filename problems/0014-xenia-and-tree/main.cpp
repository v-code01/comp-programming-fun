#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n m"; next n-1 lines = edges "u v" (1-based); next m lines =
// "t v" -- t=1 paints v red, t=2 prints edges from v to the nearest red vertex.
// stdout: one line per type-2 query.
//
// n, m up to 1e5. scanf per token is slow enough to matter at 2e5 tokens plus the
// edges, so pull raw bytes through one buffer, and stage all output in one string
// flushed once. vertex 1 is already red before the first query -- build() folds
// that in, so a leading "2 1" prints 0.
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

// every integer in this problem is a positive vertex index or a type in {1,2}.
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
    const int m = rd();

    std::vector<std::pair<int, int>> edges;
    edges.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (int i = 0; i < n - 1; ++i) {
        const int u = rd();
        const int v = rd();
        edges.emplace_back(u - 1, v - 1);  // to 0-based.
    }

    p0014::NearestRed tree;
    tree.build(n, edges);

    std::string out;
    out.reserve(static_cast<std::size_t>(m) * 3);  // rough, dodges regrowth.
    char tmp[16];

    for (int i = 0; i < m; ++i) {
        const int t = rd();
        const int v = rd() - 1;  // to 0-based.
        if (t == 1) {
            tree.paint(v);
        } else {
            const int len = std::snprintf(tmp, sizeof(tmp), "%d\n", tree.query(v));
            out.append(tmp, static_cast<std::size_t>(len));
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
