#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: "n m", then n integers, then m ops --
//   1 l r    -> print sum a[l..r]
//   2 l r x  -> a[i] %= x for i in [l, r]
//   3 k x    -> a[k] = x
// indices 1-based. one line of output per type-1 op.
//
// fast IO. m up to 1e5 with big integers on both ends -- scanf/printf per token
// would dominate the run, so we slurp stdin through a byte buffer and stage all
// output in one string flushed once. sums reach ~1e14, so int64 throughout.
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

// read one non-negative integer; every value in this problem is >= 1.
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
    int m = static_cast<int>(rd());

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    for (auto& v : a) v = rd();

    p0009::Beats tree;
    tree.build(a);

    std::string out;
    out.reserve(static_cast<std::size_t>(m) * 8);  // rough, avoids regrowth.

    for (int q = 0; q < m; ++q) {
        int type = static_cast<int>(rd());
        if (type == 1) {
            int l = static_cast<int>(rd());
            int r = static_cast<int>(rd());
            out += std::to_string(tree.range_sum(l - 1, r - 1));
            out += '\n';
        } else if (type == 2) {
            int l = static_cast<int>(rd());
            int r = static_cast<int>(rd());
            std::int64_t x = rd();
            tree.range_mod(l - 1, r - 1, x);
        } else {
            int k = static_cast<int>(rd());
            std::int64_t x = rd();
            tree.point_set(k - 1, x);
        }
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
