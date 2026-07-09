#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "solution.hpp"

// stdin: line 1 = "n q"; then q ops --
//   1 u v   -> add edge u--v
//   2 u v   -> remove edge u--v (present by contract)
//   3 u v   -> print 1 if u and v are connected right now, else 0
// vertices 1-based; one output line per type-3 op.
//
// n, q up to 2e5. all ops are read before a single answer is produced -- the whole
// method is offline. scanf per token drags at this volume, so pull raw bytes
// through one buffer and stage all output in one string flushed once.
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

// unsigned reader -- every value here (type, u, v) is a positive int.
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
    int q = rd();

    std::vector<p0039::Op> ops(static_cast<std::size_t>(q));
    for (int i = 0; i < q; ++i) {
        int type = rd();
        int u = rd();
        int v = rd();
        ops[static_cast<std::size_t>(i)] = {type, u, v};
    }

    std::vector<int> ans = p0039::solve(n, ops);

    std::string out;
    out.reserve(ans.size() * 2);  // one digit plus newline each.
    for (int a : ans) {
        out += (a ? '1' : '0');
        out += '\n';
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
