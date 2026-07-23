#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin:  line 1 = k. line 2 = 2^k values of a. line 3 = 2^k values of b.
// stdout: 2^k values of c, space-separated, where
//   c[S] = sum over T subset of S of a[T]*b[S\T]   (mod 998244353).
//
// at k=20 there are ~1e6 values per array, so scanf/printf per token is the
// bottleneck, not the O(2^k k^2) transform. slurp the input in one pass and hand
// back the output through a single buffered write.
namespace {

std::vector<char> slurp() {
    std::vector<char> buf;
    buf.reserve(1 << 22);
    char chunk[1 << 16];
    std::size_t got;
    while ((got = std::fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
        buf.insert(buf.end(), chunk, chunk + got);
    }
    return buf;
}

}  // namespace

int main() {
    const std::vector<char> buf = slurp();
    const char* p = buf.data();
    const char* const end = p + buf.size();

    // values are non-negative residues < p < 2^30, so an unsigned parse is enough.
    auto next_uint = [&]() -> std::uint32_t {
        while (p < end && (*p < '0' || *p > '9')) ++p;  // skip separators.
        std::uint32_t v = 0;
        while (p < end && *p >= '0' && *p <= '9') {
            v = v * 10u + static_cast<std::uint32_t>(*p++ - '0');
        }
        return v;
    };

    const int k = static_cast<int>(next_uint());
    const std::size_t size = static_cast<std::size_t>(1) << k;

    std::vector<std::uint32_t> a(size), b(size);
    for (std::size_t i = 0; i < size; ++i) a[i] = next_uint();
    for (std::size_t i = 0; i < size; ++i) b[i] = next_uint();

    const std::vector<std::uint32_t> c = p0056::solve(a, b, k);

    // buffered decimal writer. reserve ~11 bytes per value (10 digits + a space)
    // so the whole answer is one contiguous write, no per-number printf.
    std::vector<char> out;
    out.reserve(size * 11 + 1);
    char tmp[12];
    for (std::size_t i = 0; i < size; ++i) {
        std::uint32_t v = c[i];
        int n = 0;
        if (v == 0) {
            tmp[n++] = '0';
        } else {
            while (v > 0) {
                tmp[n++] = static_cast<char>('0' + v % 10u);
                v /= 10u;
            }
        }
        while (n > 0) out.push_back(tmp[--n]);  // digits came out reversed.
        out.push_back(i + 1 == size ? '\n' : ' ');
    }
    std::fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
