#pragma once

#include <cstdint>
#include <cstdio>

#include "common/rng.hpp"

namespace cpf {

// one tally for the whole binary. we count every break, not just the first --
// seeing all ten failures at once beats fixing one and re-running ten times.
struct TestState {
    int failed = 0;
    int total = 0;
};

inline TestState& state() {
    static TestState s;
    return s;
}

// print the tally and hand the shell an exit code it can gate on. zero is green.
inline int report() {
    TestState& s = state();
    std::printf("%d/%d checks passed\n", s.total - s.failed, s.total);
    return s.failed == 0 ? 0 : 1;
}

// the reference is dumb and obviously correct. the solution is fast and has to
// agree with it everywhere. feed both the same random input, diff the output --
// disagree on one seed and that seed is your repro, printed for you.
//
// gen(Rng&) -> Input
// sol(const Input&) -> Output
// ref(const Input&) -> Output   (Output must have operator==)
template <typename Gen, typename Sol, typename Ref>
bool differential(int cases, std::uint64_t base_seed, Gen gen, Sol sol,
                  Ref ref) {
    for (int i = 0; i < cases; ++i) {
        std::uint64_t seed = base_seed + static_cast<std::uint64_t>(i);
        Rng rng(seed);
        auto input = gen(rng);
        auto got = sol(input);
        auto want = ref(input);
        ++state().total;
        if (!(got == want)) {
            ++state().failed;
            std::printf("  DIFF at seed %llu\n",
                        static_cast<unsigned long long>(seed));
            return false;  // stop at the first break -- the seed replays it.
        }
    }
    return true;
}

}  // namespace cpf

// plain equality check with a located failure line. no expression rewriting,
// no template soup -- when it fails you want the file, the line, and the two
// names, nothing more.
#define CPF_CHECK(cond)                                                       \
    do {                                                                      \
        ++::cpf::state().total;                                               \
        if (!(cond)) {                                                        \
            ++::cpf::state().failed;                                          \
            std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);     \
        }                                                                     \
    } while (0)

#define CPF_EQ(a, b)                                                          \
    do {                                                                      \
        ++::cpf::state().total;                                               \
        if (!((a) == (b))) {                                                  \
            ++::cpf::state().failed;                                          \
            std::printf("  FAIL %s:%d  %s == %s\n", __FILE__, __LINE__, #a,   \
                        #b);                                                  \
        }                                                                     \
    } while (0)
