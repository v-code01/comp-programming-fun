#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n q"; line 2 = n values a_i; then q ops --
//   1 u v    -> link u and v (different trees)
//   2 u v    -> cut the edge u--v (it exists)
//   3 u v x  -> add x to every vertex on the path u..v
//   4 u v    -> print the sum over the path u..v
// vertices 1-based; one output line per type-4 op.
//
// n, q up to 1e5. values and x are signed and reach 1e9, path sums reach ~2.5e18
// -- int64 everywhere. scanf per token drags at this volume, so pull raw bytes
// through one buffer and stage all output in a single string flushed once.
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

// signed reader -- a_i and x can be negative, so a leading '-' has to survive.
inline std::int64_t rd() {
    int c = gc();
    while (c != -1 && c != '-' && (c < '0' || c > '9')) c = gc();
    bool neg = false;
    if (c == '-') {
        neg = true;
        c = gc();
    }
    std::int64_t x = 0;
    while (c >= '0' && c <= '9') {
        x = x * 10 + (c - '0');
        c = gc();
    }
    return neg ? -x : x;
}

}  // namespace

int main() {
    int n = static_cast<int>(rd());
    int q = static_cast<int>(rd());

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    for (auto& v : a) v = rd();

    p0035::LinkCutTree lct;
    lct.build(a);

    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 4);  // rough, dodges regrowth.

    for (int i = 0; i < q; ++i) {
        int type = static_cast<int>(rd());
        int u = static_cast<int>(rd());
        int v = static_cast<int>(rd());
        if (type == 1) {
            lct.link(u, v);
        } else if (type == 2) {
            lct.cut(u, v);
        } else if (type == 3) {
            std::int64_t x = rd();
            lct.pathAdd(u, v, x);
        } else {
            out += std::to_string(lct.pathSum(u, v));
            out += '\n';
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
