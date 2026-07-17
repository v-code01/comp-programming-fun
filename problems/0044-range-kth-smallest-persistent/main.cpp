#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n q"; line 2 = n values (|a_i| <= 1e9, so they can be
// negative); then q lines "l r k" (1-based positions, 1-based k). stdout: one
// line per query -- the k-th smallest value in a[l..r].
//
// n, q up to 1e5. scanf per token drags at that many tokens, so pull raw bytes
// through one buffer and stage all output in one string flushed once. values are
// signed, so the reader eats a leading '-'.
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

// signed base-10. skips to the first digit or '-', then folds the number.
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
    const int n = static_cast<int>(rd());
    const int q = static_cast<int>(rd());

    std::vector<std::int64_t> a(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = rd();

    p0044::PersistentKth tree;
    tree.build(a);

    std::string out;
    out.reserve(static_cast<std::size_t>(q) * 12);  // rough, dodges regrowth.
    char tmp[24];

    for (int i = 0; i < q; ++i) {
        const int l = static_cast<int>(rd());
        const int r = static_cast<int>(rd());
        const int k = static_cast<int>(rd());
        const int len = std::snprintf(tmp, sizeof(tmp), "%lld\n",
                                      static_cast<long long>(tree.query(l, r, k)));
        out.append(tmp, static_cast<std::size_t>(len));
    }

    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
