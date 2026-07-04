#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>

#include "common/rng.hpp"
#include "solution.hpp"

// the problem charges you for one suffix-array build over w = s + '#' + t, so
// |w| = 2e5 + 1 at the ceiling. two shapes:
//   all-'a' -- the doubling's worst case. every suffix ties for log n rounds, so
//   no round breaks early and the sort runs its full depth. answer is |s|.
//   random 2-letter -- a busy, near-maximal-entropy sort with short lcp runs.
// wall time is what matters; one build is tens of milliseconds, so time it once
// each rather than looping.
namespace {

double time_once(const std::string& s, const std::string& t, int& answer) {
    const auto t0 = std::chrono::steady_clock::now();
    answer = p0028::longest_common_substring(s, t);
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

}  // namespace

int main() {
    constexpr int kN = 100000;

    // all-'a': the doubling never gets to stop early.
    const std::string sa(kN, 'a');
    const std::string ta(kN, 'a');

    // random over {a,b}: dense sort, short common runs.
    cpf::Rng rng(2028);
    std::string sr(kN, 'a'), tr(kN, 'a');
    for (auto& c : sr) c = static_cast<char>('a' + rng.in_range(0, 1));
    for (auto& c : tr) c = static_cast<char>('a' + rng.in_range(0, 1));

    int ans_a = 0, ans_r = 0;
    const double ms_a = time_once(sa, ta, ans_a);
    const double ms_r = time_once(sr, tr, ans_r);

    std::printf("longest common substring  (|w| = %d)\n", 2 * kN + 1);
    std::printf("  all-'a'     %10.2f ms   answer=%d\n", ms_a, ans_a);
    std::printf("  random ab   %10.2f ms   answer=%d\n", ms_r, ans_r);
    return 0;
}
