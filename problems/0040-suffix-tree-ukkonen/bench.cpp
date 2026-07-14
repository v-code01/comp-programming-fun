#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

namespace {

// a build is one left-to-right pass -- there's no inner query to hammer, so time
// the pass itself. repeat it a handful of times and take the mean; one shot at
// 2e5 is a fraction of a millisecond and the clock noise would swamp it.
// keep() pins the answer so -O2 can't delete the build it came from.
void report(const char* name, const std::string& s) {
    constexpr int kRuns = 20;

    auto ans = p0040::suffix_tree_stats(s);  // warm: page faults, cold caches.
    cpf::keep(ans);

    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < kRuns; ++i) {
        ans = p0040::suffix_tree_stats(s);
        cpf::keep(ans);
    }
    const auto t1 = std::chrono::steady_clock::now();

    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / kRuns;
    const double mcps = static_cast<double>(s.size()) / ms / 1000.0;
    std::printf("%-10s |s|=%6zu  %7.3f ms  %7.1f M chars/s   distinct=%-12lld "
                "repeat=%lld\n",
                name, s.size(), ms, mcps,
                static_cast<long long>(ans.first),
                static_cast<long long>(ans.second));
}

}  // namespace

int main() {
    constexpr std::size_t kN = 200'000;  // the stated ceiling.

    // all-'a'. the tree is a chain of 2e5 internal nodes -- maximum depth,
    // maximum splitting, every phase walks. the deepest internal node is at
    // n-1, which is also the worst case for the dfs that reads answer (2).
    report("mono", std::string(kN, 'a'));

    // period 2. distinct substrings stay linear while the structure stays
    // maximally repetitive -- constant edge splits, constant link hops.
    {
        std::string s(kN, 'a');
        for (std::size_t i = 0; i < kN; ++i) {
            s[i] = static_cast<char>('a' + (i & 1));
        }
        report("period2", s);
    }

    // period 17. long enough that the repeats aren't trivial, short enough that
    // the active point still climbs.
    {
        std::string s(kN, 'a');
        for (std::size_t i = 0; i < kN; ++i) {
            s[i] = static_cast<char>('a' + (i % 17));
        }
        report("period17", s);
    }

    cpf::Rng rng(2025);

    // 2-letter random. repeats dense, tree deep-ish, no periodic shortcut.
    {
        std::string s(kN, 'a');
        for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, 1));
        report("rand2", s);
    }

    // 26-letter random. ~2e10 distinct substrings, all held in ~4e5 nodes --
    // the compression the whole structure is for, in one number.
    {
        std::string s(kN, 'a');
        for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, 25));
        report("rand26", s);
    }

    return 0;
}
