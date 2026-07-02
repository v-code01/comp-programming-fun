#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the two costs the problem actually charges you:
//   build -- one SAM over |s| = 1e6, the O(|s| * A) construction.
//   query -- a batch of queries whose lengths sum to ~1e6, the O(sum|x|) walk.
// pre-generate everything so the timed regions measure the automaton, not the
// rng. queries are real substrings of s, so they hit the recording+dedupe path
// instead of bailing on the first mismatch -- the honest worst case.
int main() {
    cpf::Rng rng(2026);

    // a 3-letter alphabet keeps the automaton dense: near-maximal states and a
    // busy suffix-link tree, which is where build and query both do their work.
    constexpr int kN = 1'000'000;
    std::string s(static_cast<std::size_t>(kN), 'a');
    for (auto& c : s) c = static_cast<char>('a' + rng.in_range(0, 2));

    // ---- build ----
    auto b0 = std::chrono::steady_clock::now();
    p0012::SuffixAutomaton sam;
    sam.build(s);
    auto b1 = std::chrono::steady_clock::now();
    const double build_ms =
        std::chrono::duration<double, std::milli>(b1 - b0).count();
    cpf::keep(sam.states());

    // ---- queries: random substrings of s, total length ~1e6 ----
    std::vector<std::string> queries;
    std::int64_t total_chars = 0;
    while (total_chars < kN) {
        const std::int64_t L = rng.in_range(1, 4000);
        const std::int64_t start =
            rng.in_range(0, kN - L < 0 ? 0 : kN - L);
        queries.emplace_back(s.substr(static_cast<std::size_t>(start),
                                      static_cast<std::size_t>(L)));
        total_chars += L;
    }

    auto q0 = std::chrono::steady_clock::now();
    std::int64_t sink = 0;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        sink += sam.count_cyclic(queries[i], static_cast<int>(i) + 1);
    }
    auto q1 = std::chrono::steady_clock::now();
    cpf::keep(sink);
    const double query_ms =
        std::chrono::duration<double, std::milli>(q1 - q0).count();

    std::printf("suffix automaton (alpha=3)\n");
    std::printf("  build   |s|=%d  states=%d  %10.2f ms  %8.1f Mchar/s\n", kN,
                sam.states(), build_ms, kN / build_ms / 1e3);
    std::printf("  query   %zu queries  sum|x|=%lld  %10.2f ms  %8.1f Mchar/s\n",
                queries.size(), static_cast<long long>(total_chars), query_ms,
                static_cast<double>(total_chars) / query_ms / 1e3);
    std::printf("  (sink=%lld -- keeps the optimizer honest)\n",
                static_cast<long long>(sink));
    return 0;
}
