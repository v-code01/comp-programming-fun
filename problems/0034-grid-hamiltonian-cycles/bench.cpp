#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/bench.hpp"
#include "solution.hpp"

namespace {

std::vector<std::string> full(int n, int m) {
    return std::vector<std::string>(static_cast<std::size_t>(n),
                                    std::string(static_cast<std::size_t>(m), '.'));
}

}  // namespace

// the worst cases the constraints allow: full 12x12 and 10x12, no obstacles, so
// the profile is at its widest and every cell is a live transition. the number
// printed is one full solve -- the plug DP sweeping 144 (and 120) cells over the
// whole reachable bracket-profile set.
int main() {
    const std::vector<std::string> g12 = full(12, 12);
    const std::vector<std::string> g10 = full(10, 12);

    cpf::bench("ham-cycles full 12x12", 5,
               [&] { return p0034::count_hamiltonian_cycles(g12); });
    cpf::bench("ham-cycles full 10x12", 5,
               [&] { return p0034::count_hamiltonian_cycles(g10); });

    std::printf(
        "(ns/op above is one full solve; /1e6 = ms. answers are mod 1e9+7.)\n");
    std::printf("12x12 count mod 1e9+7 = %lld\n",
                static_cast<long long>(p0034::count_hamiltonian_cycles(g12)));
    std::printf("10x12 count mod 1e9+7 = %lld\n",
                static_cast<long long>(p0034::count_hamiltonian_cycles(g10)));
    return 0;
}
