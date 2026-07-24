#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n q"; then q ops, processed strictly in order -- ONLINE --
//   1 u v w  -> link: add edge {u, v} of weight w (u, v in different trees)
//   2 u v    -> cut: remove the edge {u, v} (it exists)
//   3 u      -> print the diameter of the tree containing u (0 for a lone vertex)
// vertices 1-based; one output line per type-3 op.
//
// n up to 1e5, q up to 2e5, w up to 1e9, a diameter up to ~1e14 -- int64 throughout.
// scanf per token drags at this volume, so pull raw bytes through one buffer and
// stage all output in one string flushed once.
namespace {

constexpr int kBuf = 1 << 16;
char g_in[kBuf];
int g_pos = 0;
int g_len = 0;

inline int gc() {
    if (g_pos == g_len) {
        g_len = static_cast<int>(std::fread(g_in, 1, kBuf, stdin));
        g_pos = 0;
        if (g_len == 0) return -1;  // EOF.
    }
    return g_in[g_pos++];
}

// every value here (n, q, type, u, v, w) is a positive integer, w reaches 1e9.
inline std::int64_t rd() {
    int c = gc();
    while (c != -1 && (c < '0' || c > '9')) c = gc();
    std::int64_t x = 0;
    while (c >= '0' && c <= '9') {
        x = x * 10 + (c - '0');
        c = gc();
    }
    return x;
}

}  // namespace

int main() {
    int n = static_cast<int>(rd());
    int q = static_cast<int>(rd());

    p0058::DynamicDiameter dd(n);

    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 4);  // rough, dodges regrowth.

    for (int i = 0; i < q; ++i) {
        int type = static_cast<int>(rd());
        if (type == 1) {
            int u = static_cast<int>(rd());
            int v = static_cast<int>(rd());
            std::int64_t w = rd();
            dd.link(u, v, w);
        } else if (type == 2) {
            int u = static_cast<int>(rd());
            int v = static_cast<int>(rd());
            dd.cut(u, v);
        } else {
            int u = static_cast<int>(rd());
            out += std::to_string(dd.diameter(u));
            out += '\n';
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
