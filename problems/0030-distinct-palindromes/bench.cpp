#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

// one full build, timed once. the eertree is a single left-to-right pass, so a
// repeated micro-bench loop would just re-measure the same 1e6-char pass -- time
// it once and report the wall clock. keep() pins the result so -O2 can't delete
// the build.
double time_ms(const std::string& s, std::int64_t& out) {
    auto t0 = std::chrono::steady_clock::now();
    out = p0030::count_distinct(s);
    auto t1 = std::chrono::steady_clock::now();
    cpf::keep(out);
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void report(const char* name, const std::string& s) {
    std::int64_t ans = 0;
    const double ms = time_ms(s, ans);
    const double mcps = static_cast<double>(s.size()) / ms / 1000.0;
    std::printf("%-8s |s|=%zu  %8.2f ms  %10lld distinct  %7.1f M chars/s\n",
                name, s.size(), ms, static_cast<long long>(ans), mcps);
}

}  // namespace

int main() {
    constexpr std::size_t kN = 1'000'000;

    // worst case for the output: all one letter. every length is a brand-new
    // palindrome, so the tree grows a node on every single character and the walk
    // stays hot. the answer is |s| itself -- the linear lower bound, made concrete.
    report("mono", std::string(kN, 'a'));

    cpf::Rng rng(2025);

    // 2-letter random -- palindromes dense, plenty of link climbing.
    {
        std::string s(kN, 'a');
        for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, 1));
        report("rand2", s);
    }

    // full 26-letter random -- palindromes sparse, short walks, the easy end.
    {
        std::string s(kN, 'a');
        for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, 25));
        report("rand26", s);
    }

    return 0;
}
