#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// the cost the problem actually charges: two FFT convolutions over a transform
// sized to |T|+|P|. bench the stated worst case, |T|=2e5 and |P|=1e5, so the
// transform is 2^19 = 524288. a 3-letter alphabet plus wildcards keeps matches
// plentiful -- the position collection loop does real work instead of skipping.
// build the strings first so the timed region is the matcher, not the rng.
int main() {
    cpf::Rng rng(2026);

    constexpr int kN = 200000;
    constexpr int kM = 100000;

    std::string text(static_cast<std::size_t>(kN), 'a');
    for (auto& c : text)
        c = static_cast<char>('a' + rng.in_range(0, 2));

    std::string pat(static_cast<std::size_t>(kM), 'a');
    for (auto& c : pat) {
        if (rng.in_range(0, 99) < 20)
            c = '?';  // ~20% wildcards.
        else
            c = static_cast<char>('a' + rng.in_range(0, 2));
    }

    // warm once (cold i-cache, first page faults, roots-table build), then time.
    volatile std::size_t sink = 0;
    {
        const std::vector<int> hits = p0031::matches(text, pat);
        sink += hits.size();
    }

    constexpr int kIters = 20;
    auto t0 = std::chrono::steady_clock::now();
    std::size_t last = 0;
    for (int it = 0; it < kIters; ++it) {
        const std::vector<int> hits = p0031::matches(text, pat);
        last = hits.size();
        sink += last;
    }
    auto t1 = std::chrono::steady_clock::now();
    cpf::keep(sink);

    const double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / kIters;
    const double nplusm = static_cast<double>(kN) + static_cast<double>(kM);
    std::printf("wildcard FFT match\n");
    std::printf("  |T|=%d  |P|=%d  transform=524288\n", kN, kM);
    std::printf("  %8.2f ms/run   %8.1f Mchar/s   matches=%zu\n", ms,
                nplusm / ms / 1e3, last);
    std::printf("  (sink=%zu -- keeps the optimizer honest)\n",
                static_cast<std::size_t>(sink));
    return 0;
}
