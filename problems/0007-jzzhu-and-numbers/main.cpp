#include <cstdint>
#include <cstdio>
#include <vector>

#include "solution.hpp"

// stdin: n, then n integers in [0, 2^20). stdout: the count of non-empty groups
// whose bitwise AND is 0, mod 1e9+7.
//
// n reaches 1e6, so scanf per token is the bottleneck, not the SOS. slurp the
// whole stream once and hand-parse unsigned ints -- input drops to a single pass.
namespace {

// fixed 20-bit domain: a_i < 2^20. matches the SOS array the solver builds.
constexpr int kBits = 20;

// read the entire stream into one buffer. one syscall path, no per-token cost.
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
    std::vector<char> buf = slurp();
    const char* p = buf.data();
    const char* end = p + buf.size();

    auto next_uint = [&]() -> std::uint32_t {
        while (p < end && (*p < '0' || *p > '9')) ++p;  // skip any separators.
        std::uint32_t v = 0;
        while (p < end && *p >= '0' && *p <= '9') v = v * 10 + (*p++ - '0');
        return v;
    };

    std::uint32_t n = next_uint();
    std::vector<std::uint32_t> a(n);
    for (std::uint32_t i = 0; i < n; ++i) a[i] = next_uint();

    std::printf("%lld\n", static_cast<long long>(p0007::solve(a, kBits)));
    return 0;
}
