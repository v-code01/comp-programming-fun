#pragma once

#include <cstdint>
#include <random>
#include <vector>

namespace cpf {

// seeded on purpose. a differential failure hands you back the seed, and the
// seed replays the exact case -- randomness you can't reproduce is randomness
// you can't debug.
class Rng {
public:
    explicit Rng(std::uint64_t seed) : gen_(seed) {}

    // both ends inclusive. contest bounds are usually written [lo, hi], so
    // match them and skip the off-by-one arithmetic at every call site.
    std::int64_t in_range(std::int64_t lo, std::int64_t hi) {
        std::uniform_int_distribution<std::int64_t> d(lo, hi);
        return d(gen_);
    }

    std::vector<std::int64_t> vec(std::size_t n, std::int64_t lo,
                                  std::int64_t hi) {
        std::vector<std::int64_t> v(n);
        for (auto& x : v) x = in_range(lo, hi);
        return v;
    }

private:
    // mersenne twister -- overkill for stream quality, but the 64-bit variant
    // is fast enough that generation never shows up in the profile.
    std::mt19937_64 gen_;
};

}  // namespace cpf
