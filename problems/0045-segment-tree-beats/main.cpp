#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: "n q", then n integers, then q ops --
//   1 l r x  -> a[i] = min(a[i], x) for i in [l, r]
//   2 l r x  -> a[i] += x           for i in [l, r]   (x may be negative)
//   3 l r    -> print max of a[l..r]
//   4 l r    -> print sum of a[l..r]
// indices 1-based. one output line per type-3 and type-4 op.
//
// fast IO. q up to 1e5 with signed values on both ends -- scanf/printf per token
// would dominate, so we slurp stdin through a byte buffer and stage all output
// in one string flushed once. values and adds are signed; sums reach ~1e18 under
// the statement's magnitude cap, so int64 throughout.
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

// read one signed integer -- values, adds, and chmin caps can all be negative.
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

    p0045::Beats tree;
    tree.build(a);

    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 8);  // rough, avoids regrowth.

    for (int i = 0; i < q; ++i) {
        int type = static_cast<int>(rd());
        int l = static_cast<int>(rd());
        int r = static_cast<int>(rd());
        if (type == 1) {
            std::int64_t x = rd();
            tree.chmin(l - 1, r - 1, x);
        } else if (type == 2) {
            std::int64_t x = rd();
            tree.add(l - 1, r - 1, x);
        } else if (type == 3) {
            out += std::to_string(tree.query_max(l - 1, r - 1));
            out += '\n';
        } else {
            out += std::to_string(tree.query_sum(l - 1, r - 1));
            out += '\n';
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
