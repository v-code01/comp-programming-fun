#include <cstdint>
#include <string>

#include "common/bench.hpp"
#include "common/rng.hpp"
#include "solution.hpp"

// worst case is n = 500 -- the constraint ceiling. the inner k-scan runs the
// full interval every time regardless of the letters, so the O(n^3) count is
// fixed at ~2e7 dp writes each doing up to n work; content only steers which
// branch of the min wins. two shapes to be sure:
//   "abab…" -- two letters, so nearly every k hits s[k]==s[l] and the merge
//   term fires on every step. the busy case.
//   26 random letters -- the merge rarely fires, mostly the +1 fallback.
// same asymptotic work, and the numbers say the content barely moves the wall
// clock. iters are low because one solve is tens of milliseconds -- that's the
// point of the benchmark.
int main() {
    constexpr int kN = 500;

    std::string alt(kN, 'a');
    for (int i = 0; i < kN; ++i) alt[static_cast<std::size_t>(i)] =
        static_cast<char>('a' + (i & 1));

    cpf::Rng rng(2025);
    std::string rnd(kN, 'a');
    for (auto& c : rnd) c = static_cast<char>('a' + rng.in_range(0, 25));

    cpf::bench("clear n=500 abab", 40, [&] { return p0006::clear(alt); });
    cpf::bench("clear n=500 rand26", 40, [&] { return p0006::clear(rnd); });
    return 0;
}
