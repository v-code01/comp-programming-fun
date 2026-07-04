#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the two costs the problem charges: build the automaton over sum|P| = 1e5, then
// run the text of |T| = 1e5 through it once. pre-generate everything so the timed
// regions measure the automaton, not the rng. patterns are real substrings of T,
// so the text feed lands on live terminals -- the busy path, the honest case.
int main() {
    cpf::Rng rng(2029);

    // a 3-letter alphabet keeps the trie dense: near-maximal branching and a busy
    // fail tree, where both build and feed do their work.
    constexpr int kT = 100'000;
    constexpr int kSumP = 100'000;

    std::string text(static_cast<std::size_t>(kT), 'a');
    for (auto& c : text) c = static_cast<char>('a' + rng.in_range(0, 2));

    // patterns: real substrings of T, lengths summing to ~sum|P|.
    std::vector<std::string> patterns;
    std::int64_t total_chars = 0;
    while (total_chars < kSumP) {
        const std::int64_t L = rng.in_range(1, 40);
        const std::int64_t start = rng.in_range(0, kT - L < 0 ? 0 : kT - L);
        patterns.emplace_back(text.substr(static_cast<std::size_t>(start),
                                          static_cast<std::size_t>(L)));
        total_chars += L;
    }

    // ---- build ----
    auto b0 = std::chrono::steady_clock::now();
    p0029::AhoCorasick ac;
    ac.build(patterns);
    auto b1 = std::chrono::steady_clock::now();
    const double build_ms =
        std::chrono::duration<double, std::milli>(b1 - b0).count();
    cpf::keep(ac.states());

    // ---- feed: the O(|T|) text run, one goto read and one add per char ----
    auto f0 = std::chrono::steady_clock::now();
    const std::int64_t answer = ac.count(text);
    auto f1 = std::chrono::steady_clock::now();
    cpf::keep(answer);
    const double feed_ms =
        std::chrono::duration<double, std::milli>(f1 - f0).count();

    std::printf("aho-corasick (alpha=3)\n");
    std::printf("  build   sum|P|=%lld  patterns=%zu  states=%d  %8.2f ms  "
                "%8.1f Mchar/s\n",
                static_cast<long long>(total_chars), patterns.size(),
                ac.states(), build_ms,
                static_cast<double>(total_chars) / build_ms / 1e3);
    std::printf("  feed    |T|=%d  %8.2f ms  %8.1f Mchar/s\n", kT, feed_ms,
                kT / feed_ms / 1e3);
    std::printf("  total occurrences=%lld\n", static_cast<long long>(answer));
    return 0;
}
