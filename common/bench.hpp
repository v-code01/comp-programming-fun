#pragma once

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>

namespace cpf {

// stop the optimizer from deleting work whose result we throw away. under -O2
// the compiler sees a loop that computes nothing observable and erases it --
// you'd "benchmark" an empty loop and report 0 ns. this pins the value to a
// register or memory so the work has to happen.
template <typename T>
inline void keep(const T& value) {
    asm volatile("" : : "r,m"(value) : "memory");
}

// warm first, then time. the first runs pay for page faults, cold i-cache, and
// an untrained branch predictor -- costs you never pay again in a hot loop.
// timing them in slanders the code.
template <typename F>
void bench(const std::string& name, std::uint64_t iters, F&& f) {
    std::uint64_t warm = iters / 10 + 1;  // 10% thrown away, minimum one.
    for (std::uint64_t i = 0; i < warm; ++i) keep(f());

    auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) keep(f());
    auto t1 = std::chrono::steady_clock::now();

    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    double per = ns / static_cast<double>(iters);
    std::printf("%-28s %12.2f ns/op  %14.0f ops/s\n", name.c_str(), per,
                1e9 / per);
}

}  // namespace cpf
