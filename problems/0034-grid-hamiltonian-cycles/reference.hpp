#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace p0034 {

// the dumb, obviously-correct oracle -- a depth-first walk that shares NOT ONE
// idea with the plug DP. no profiles, no brackets, no connectivity bookkeeping.
// it just walks the loop.
//
// pick the first free cell as the start and fix it. from there, step to any
// adjacent free cell you haven't stood on yet. when every free cell has been
// stood on, ask one thing: is the cell you're on adjacent to the start? if yes,
// the walk closed into a loop that touched every free cell exactly once -- a
// hamiltonian cycle.
//
// each undirected cycle gets walked twice, once each direction (start -> a ->
// … -> b -> start and its mirror). so count directed closed walks and divide by
// two. no reflection to reason about, no first-move trick -- just halve at the
// end.
//
// fewer than four free cells can't hold a grid cycle: the grid is bipartite so a
// cycle needs an even length, and the shortest even simple cycle here is the 2x2
// square of four cells. bail to 0 below four -- and that bail also kills the
// degenerate "walk one edge and back" that a raw adjacency check would miscount.
//
// exponential in the free-cell count -- that's the point. it's the ground truth
// for small grids, nothing more. keep the free-cell count at ~16 or under.
class ReferenceSolver {
public:
    std::int64_t solve(const std::vector<std::string>& grid) {
        n_ = static_cast<int>(grid.size());
        if (n_ == 0) return 0;
        m_ = static_cast<int>(grid[0].size());

        // number every free cell and remember where it sits.
        id_.assign(static_cast<std::size_t>(n_),
                   std::vector<int>(static_cast<std::size_t>(m_), -1));
        cells_.clear();
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j)
                if (grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] ==
                    '.') {
                    id_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                        static_cast<int>(cells_.size());
                    cells_.push_back({i, j});
                }

        free_ = static_cast<int>(cells_.size());
        if (free_ < 4) return 0;  // no grid cycle below four cells.

        // adjacency list over free cells, 4-directional.
        adj_.assign(static_cast<std::size_t>(free_), {});
        for (int c = 0; c < free_; ++c) {
            const auto [i, j] = cells_[static_cast<std::size_t>(c)];
            static constexpr std::array<int, 4> di{-1, 1, 0, 0};
            static constexpr std::array<int, 4> dj{0, 0, -1, 1};
            for (int d = 0; d < 4; ++d) {
                const int ni = i + di[static_cast<std::size_t>(d)];
                const int nj = j + dj[static_cast<std::size_t>(d)];
                if (ni < 0 || ni >= n_ || nj < 0 || nj >= m_) continue;
                const int nid =
                    id_[static_cast<std::size_t>(ni)][static_cast<std::size_t>(nj)];
                if (nid >= 0)
                    adj_[static_cast<std::size_t>(c)].push_back(nid);
            }
        }

        // any free cell with fewer than two free neighbors can't have degree two
        // in a cycle -- a dead end. no need to search, the answer is 0.
        for (int c = 0; c < free_; ++c)
            if (adj_[static_cast<std::size_t>(c)].size() < 2) return 0;

        vis_.assign(static_cast<std::size_t>(free_), 0);
        vis_[0] = 1;            // start is cell 0, fixed.
        walks_ = 0;
        dfs(0, 1);             // at cell 0, one cell visited so far.
        return walks_ / 2;     // two directions per undirected cycle.
    }

private:
    void dfs(int cur, int seen) {
        if (seen == free_) {
            // every cell stood on -- close it iff we can step back to start.
            for (int nx : adj_[static_cast<std::size_t>(cur)])
                if (nx == 0) ++walks_;
            return;
        }
        for (int nx : adj_[static_cast<std::size_t>(cur)]) {
            if (vis_[static_cast<std::size_t>(nx)]) continue;
            vis_[static_cast<std::size_t>(nx)] = 1;
            dfs(nx, seen + 1);
            vis_[static_cast<std::size_t>(nx)] = 0;
        }
    }

    int n_ = 0, m_ = 0, free_ = 0;
    std::int64_t walks_ = 0;
    std::vector<std::pair<int, int>> cells_;
    std::vector<std::vector<int>> id_;
    std::vector<std::vector<int>> adj_;
    std::vector<char> vis_;
};

// free function wrapper -- matches the solution's call shape.
inline std::int64_t count_hamiltonian_cycles_ref(
    const std::vector<std::string>& grid) {
    ReferenceSolver s;
    return s.solve(grid);
}

}  // namespace p0034
