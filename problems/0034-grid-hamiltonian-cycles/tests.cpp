#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "common/rng.hpp"
#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

namespace {

using ll = std::int64_t;

ll fast(const std::vector<std::string>& g) {
    return p0034::count_hamiltonian_cycles(g);
}
ll brute(const std::vector<std::string>& g) {
    return p0034::count_hamiltonian_cycles_ref(g);
}

// a full n x m grid, no obstacles.
std::vector<std::string> full(int n, int m) {
    return std::vector<std::string>(static_cast<std::size_t>(n),
                                    std::string(static_cast<std::size_t>(m), '.'));
}

// a random grid. each cell free with probability free_pct/100. small dims keep
// the brute oracle's exponential walk cheap.
std::vector<std::string> random_grid(cpf::Rng& rng, int nlo, int nhi, int mlo,
                                     int mhi, int free_pct) {
    const int n = static_cast<int>(rng.in_range(nlo, nhi));
    const int m = static_cast<int>(rng.in_range(mlo, mhi));
    std::vector<std::string> g(static_cast<std::size_t>(n),
                               std::string(static_cast<std::size_t>(m), '#'));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (rng.in_range(1, 100) <= free_pct)
                g[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = '.';
    return g;
}

int count_free(const std::vector<std::string>& g) {
    int f = 0;
    for (const auto& row : g)
        for (char c : row)
            if (c == '.') ++f;
    return f;
}

}  // namespace

int main() {
    using namespace p0034;

    // ---- the OEIS cross-checks. every "expected" is the real count AND agrees
    // with the brute walk, so no number here is asserted from memory alone. ----

    // 2x2 -- the single square. 1 cycle.
    CPF_EQ(fast(full(2, 2)), 1);
    CPF_EQ(brute(full(2, 2)), 1);

    // 2x3 -- only the outer boundary loops. 1.
    CPF_EQ(fast(full(2, 3)), 1);
    CPF_EQ(brute(full(2, 3)), 1);

    // 2x4 -- again just the boundary. 1.
    CPF_EQ(fast(full(2, 4)), 1);
    CPF_EQ(brute(full(2, 4)), 1);

    // 3x3 -- nine cells, odd. the grid is bipartite, a cycle needs equal colors,
    // so no hamiltonian cycle exists. 0.
    CPF_EQ(fast(full(3, 3)), 0);
    CPF_EQ(brute(full(3, 3)), 0);

    // 4x4 -- the textbook value, A003763(2). 6.
    CPF_EQ(fast(full(4, 4)), 6);
    CPF_EQ(brute(full(4, 4)), 6);

    // a couple more full grids, solution cross-checked against the brute walk.
    CPF_EQ(fast(full(2, 12)), brute(full(2, 12)));  // long thin strip -> 1.
    CPF_EQ(fast(full(2, 12)), 1);
    CPF_EQ(fast(full(4, 6)), brute(full(4, 6)));
    CPF_EQ(fast(full(3, 4)), brute(full(3, 4)));
    CPF_EQ(fast(full(4, 4)), brute(full(4, 4)));

    // ---- edges ----

    // everything blocked -> nothing to loop. 0.
    {
        std::vector<std::string> g = {"###", "###"};
        CPF_EQ(fast(g), 0);
        CPF_EQ(brute(g), 0);
    }

    // one lone free cell -> can't have degree two. 0.
    {
        std::vector<std::string> g = {"###", "#.#", "###"};
        CPF_EQ(fast(g), 0);
        CPF_EQ(brute(g), 0);
    }

    // a 2x2 hole in the middle of a 4x4 -> the 12-cell border ring, exactly one
    // cycle.
    {
        std::vector<std::string> g = {"....", ".##.", ".##.", "...."};
        CPF_EQ(fast(g), 1);
        CPF_EQ(brute(g), 1);
    }

    // an L-shaped region. the elbow forces the walk; whatever the count, both
    // agree.
    {
        std::vector<std::string> g = {"..#", "..#", "..."};
        CPF_EQ(fast(g), brute(g));
    }

    // a wall splits the grid into two 2x2 rooms -> two loops needed, never one.
    // 0.
    {
        std::vector<std::string> g = {"..#..", "..#..", "..#.."};
        CPF_EQ(fast(g), 0);
        CPF_EQ(brute(g), 0);
    }

    // a free cell with a single free neighbor -- a spur off a room. dead end, 0.
    {
        std::vector<std::string> g = {"...", ".#.", "#.#"};  // (2,1) has one nbr.
        CPF_EQ(fast(g), 0);
        CPF_EQ(brute(g), 0);
    }

    // a plus/cross shape -- a middle cell of degree four, four arms of one cell.
    // each arm cell has a single neighbor, so no cycle. 0.
    {
        std::vector<std::string> g = {"#.#", "...", "#.#"};
        CPF_EQ(fast(g), 0);
        CPF_EQ(brute(g), 0);
    }

    // a 3x4 full grid -- has real cycles, both engines land on the same number.
    {
        CPF_EQ(fast(full(3, 4)), brute(full(3, 4)));
    }

    // ---- randomized differential vs the brute walk. thousands of small grids
    // with obstacles, so the region shape and every merge/close transition gets
    // hammered. small dims keep the oracle honest and fast. ----

    int total_cases = 0, nonzero = 0, diffs = 0, first_seed = -1;

    // sweep A -- tiny grids, medium obstacle density. dims 2..4, ~65% free.
    {
        cpf::Rng rng(0x0034A);
        const int kCases = 40000;
        for (int c = 0; c < kCases; ++c) {
            auto g = random_grid(rng, 2, 4, 2, 4, 65);
            if (count_free(g) > 16) continue;  // keep the brute walk cheap.
            const ll got = fast(g), want = brute(g);
            ++total_cases;
            if (want > 0) ++nonzero;
            if (got != want && diffs == 0) first_seed = c;
            if (got != want) ++diffs;
        }
    }

    // sweep B -- denser, so full rooms and real loops show up often. dims 2..5,
    // ~78% free, capped at 16 free cells.
    {
        cpf::Rng rng(0x0034B);
        const int kCases = 40000;
        for (int c = 0; c < kCases; ++c) {
            auto g = random_grid(rng, 2, 5, 2, 5, 78);
            if (count_free(g) > 16) continue;
            const ll got = fast(g), want = brute(g);
            ++total_cases;
            if (want > 0) ++nonzero;
            if (got != want && diffs == 0) first_seed = c;
            if (got != want) ++diffs;
        }
    }

    // sweep C -- sparse holes punched into otherwise-full grids. start full,
    // knock out a few cells -> lots of ring-with-holes shapes.
    {
        cpf::Rng rng(0x0034C);
        const int kCases = 30000;
        for (int c = 0; c < kCases; ++c) {
            const int n = static_cast<int>(rng.in_range(3, 4));
            const int m = static_cast<int>(rng.in_range(3, 5));
            auto g = full(n, m);
            const int holes = static_cast<int>(rng.in_range(0, 4));
            for (int h = 0; h < holes; ++h) {
                const int hi = static_cast<int>(rng.in_range(0, n - 1));
                const int hj = static_cast<int>(rng.in_range(0, m - 1));
                g[static_cast<std::size_t>(hi)][static_cast<std::size_t>(hj)] = '#';
            }
            if (count_free(g) > 16) continue;
            const ll got = fast(g), want = brute(g);
            ++total_cases;
            if (want > 0) ++nonzero;
            if (got != want && diffs == 0) first_seed = c;
            if (got != want) ++diffs;
        }
    }

    // sweep D -- non-square, so the transpose-to-min-dimension path runs both
    // ways. tall grids (more rows than cols) and wide grids both appear.
    {
        cpf::Rng rng(0x0034D);
        const int kCases = 20000;
        for (int c = 0; c < kCases; ++c) {
            auto g = random_grid(rng, 2, 6, 2, 3, 80);  // tall and narrow.
            if (count_free(g) > 16) continue;
            const ll got = fast(g), want = brute(g);
            ++total_cases;
            if (want > 0) ++nonzero;
            if (got != want && diffs == 0) first_seed = c;
            if (got != want) ++diffs;
        }
    }

    CPF_CHECK(diffs == 0);
    if (diffs) std::printf("  first diff in a sweep at case %d\n", first_seed);
    std::printf(
        "differential vs brute walk: %d cases, %d nonzero, %d diffs\n",
        total_cases, nonzero, diffs);

    return cpf::report();
}
