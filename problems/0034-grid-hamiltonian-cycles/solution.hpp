#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace p0034 {

constexpr std::int64_t MOD = 1000000000LL + 7;

// count hamiltonian cycles on a grid with obstacles -- closed loops that pass
// through every free cell exactly once, moving only between orthogonally
// adjacent free cells, forming ONE cycle, not several disjoint loops. answer mod
// 1e9+7. this is the classic plug DP, "formula 1" (URAL 1519).
//
// ---- the idea ----
// sweep cell by cell, row-major. draw the contour that splits processed cells
// from the rest -- a broken staircase. the loop's edges that CROSS that contour
// are the "plugs." track them, and you never need to remember the interior: the
// contour is the whole state.
//
// across the contour there are exactly m+1 plug slots. about to process (i, j):
//   slots 0 .. j-1   -- the down-edge under each cell (i, k), k < j.
//   slot  j          -- the horizontal edge between (i, j-1) and (i, j). the
//                        LEFT plug of the current cell.
//   slots j+1 .. m   -- the down-edge under each cell (i-1, k), k >= j. slot j+1
//                        is the down-edge under (i-1, j) -- the UP plug of the
//                        current cell.
// so processing (i, j) reads slot j (left) and slot j+1 (up), and writes slot j
// (the cell's new DOWN plug) and slot j+1 (its new RIGHT plug). after the last
// column, shift the whole array one slot right and drop slot m -- that's the
// step to the next row.
//
// ---- the plug encoding: brackets ----
// a plug isn't just present/absent. two plugs that are the two dangling ends of
// the same half-built path have to KNOW they're paired -- close the wrong two and
// you seal a short sub-loop. the loop is planar, so its crossings nest like
// parentheses and never tangle. that lets three states per plug carry the whole
// connectivity:
//   0  -- no plug. no edge crosses here.
//   1  -- "(", the LEFT end of a component. its partner is a matching ")" to the
//         right.
//   2  -- ")", the RIGHT end. its partner is a matching "(" to the left.
// two plugs are joined iff they're a matched bracket pair. the profile is always
// a balanced, properly-nested bracket string (with 0s as gaps). pack it two bits
// per slot into one integer and hash on it.
//
// ---- the four transition shapes at a free cell (L = left plug, U = up plug) ----
// a free cell needs degree two in the loop. left and up are fixed by the
// profile; right and down are ours to place.
//
//   1. both empty (L=0, U=0). the loop hasn't reached this cell. it must OPEN a
//      brand-new component here, using its right and down edges. down is the left
//      end, right the right end -- write "(" to slot j, ")" to slot j+1. needs a
//      free right neighbor and a free down neighbor, else this cell is a dead end.
//
//   2. exactly one present (L xor U). the cell EXTENDS that one path with a
//      single outgoing edge, right or down. the bracket value just rides along to
//      its new slot -- an open stays open, a close stays close. the other slot
//      goes empty.
//
//   3. both present, and they DON'T seal the same component (L=")" U="(", or two
//      opens, or two closes). the cell MERGES two components -- it spends both its
//      left and up edges, writes 0 to slots j and j+1, and re-pairs the brackets
//      so the surviving two ends read as one component:
//        L="(" U="(" -> the inner "(" (U's match, forward from j+1) had a partner
//                       ")"; flip that ")" to "(", it's now the merged left end.
//        L=")" U=")" -> flip L's partner "(" (backward from j) to ")", the merged
//                       right end.
//        L=")" U="(" -> the two ends already read left-then-right. remove both,
//                       flip nothing.
//
//   4. both present AND they seal the same component (L="(", U=")"). adjacency
//      forces "(" at j to match ")" at j+1, so joining them CLOSES a loop. this
//      is the whole trap.
//
// ---- WHY early-close is forbidden ----
// a close welds a component's own two ends into a finished loop. that loop can
// never grow again -- it's shut. if any other plug is still open when you close,
// or if any free cell still lies ahead, that other stuff has to become its own
// separate loop. two loops is not a hamiltonian cycle. so a close is legal in
// exactly one spot: the LAST free cell, and only when removing these two ends
// leaves the profile completely empty. that lone close is the single global
// cycle sealing shut -- accumulate the count there. everywhere else a close is a
// dead branch, dropped on the floor.
//
// ---- complexity and the lower bound ----
// O(n * m * S) where S is the number of reachable bracket profiles of width m+1.
// S is Catalan-like, ~O(3^m) -- exponential in the profile width, which is why m
// is pinned to min(n, m) <= 12 by transposing. reading the grid is Omega(n*m),
// so the input floor is polynomial, but the work is exponential in m and no
// polynomial algorithm is known: counting hamiltonian cycles is #P-hard in
// general. keeping the profile narrow is the only lever there is.
inline std::int64_t count_hamiltonian_cycles(
    const std::vector<std::string>& grid_in) {
    int n = static_cast<int>(grid_in.size());
    if (n == 0) return 0;
    int m = static_cast<int>(grid_in[0].size());
    if (m == 0) return 0;

    // run the profile over the smaller dimension -- state count is exponential in
    // it. transposing leaves the cycle count untouched (a loop is a loop under a
    // 90-degree turn), so make columns the short side.
    std::vector<std::string> grid;
    if (n < m) {
        grid.assign(static_cast<std::size_t>(m),
                    std::string(static_cast<std::size_t>(n), '.'));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                grid[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] =
                    grid_in[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        std::swap(n, m);
    } else {
        grid = grid_in;
    }

    // find the last free cell in sweep order -- the only place a close may seal
    // the global cycle. also count free cells: below four there's no cycle.
    int last_i = -1, last_j = -1, free_cnt = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] ==
                '.') {
                ++free_cnt;
                last_i = i;
                last_j = j;
            }
    if (free_cnt < 4) return 0;

    // two bits per plug, m+1 plugs. slot k lives in bits [2k, 2k+1].
    const std::uint64_t mask =
        (static_cast<std::uint64_t>(1) << (2 * (m + 1))) - 1;

    auto get_plug = [](std::uint64_t code, int k) -> int {
        return static_cast<int>((code >> (2 * k)) & 3u);
    };
    auto set_plug = [](std::uint64_t code, int k, int v) -> std::uint64_t {
        code &= ~(static_cast<std::uint64_t>(3) << (2 * k));
        code |= static_cast<std::uint64_t>(v) << (2 * k);
        return code;
    };
    // match of a "(" at slot s -- scan right, the first ")" at equal depth.
    auto match_forward = [&](std::uint64_t code, int s) -> int {
        int depth = 0;
        for (int k = s + 1; k <= m; ++k) {
            const int v = get_plug(code, k);
            if (v == 1)
                ++depth;
            else if (v == 2) {
                if (depth == 0) return k;
                --depth;
            }
        }
        return -1;  // balanced invariant -- unreachable.
    };
    // match of a ")" at slot s -- scan left, the first "(" at equal depth.
    auto match_backward = [&](std::uint64_t code, int s) -> int {
        int depth = 0;
        for (int k = s - 1; k >= 0; --k) {
            const int v = get_plug(code, k);
            if (v == 2)
                ++depth;
            else if (v == 1) {
                if (depth == 0) return k;
                --depth;
            }
        }
        return -1;  // unreachable.
    };

    std::unordered_map<std::uint64_t, std::int64_t> cur;
    cur.reserve(1 << 12);
    cur.emplace(0, 1);  // empty contour, one way to have drawn nothing.

    std::int64_t answer = 0;

    auto add = [](std::unordered_map<std::uint64_t, std::int64_t>& dst,
                  std::uint64_t code, std::int64_t v) {
        auto& r = dst[code];
        r = (r + v) % MOD;
    };

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            const bool free_cell =
                grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] ==
                '.';
            const bool last_cell = (i == last_i && j == last_j);
            const bool right_free =
                (j + 1 < m) &&
                grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j + 1)] ==
                    '.';
            const bool down_free =
                (i + 1 < n) &&
                grid[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(j)] ==
                    '.';

            std::unordered_map<std::uint64_t, std::int64_t> nxt;
            nxt.reserve(cur.size() * 2);

            for (const auto& [code, cnt] : cur) {
                const int L = get_plug(code, j);
                const int U = get_plug(code, j + 1);

                if (!free_cell) {
                    // a wall. no plug may run into it -- a path can't dead-end at
                    // an obstacle. valid only when both slots are empty; it then
                    // propagates the profile untouched.
                    if (L == 0 && U == 0) add(nxt, code, cnt);
                    continue;
                }

                if (L == 0 && U == 0) {
                    // shape 1 -- open a new component with right and down.
                    if (right_free && down_free) {
                        std::uint64_t nc = set_plug(code, j, 1);   // down = "("
                        nc = set_plug(nc, j + 1, 2);               // right = ")"
                        add(nxt, nc, cnt);
                    }
                } else if (L != 0 && U == 0) {
                    // shape 2 -- extend the left path, carrying its bracket.
                    if (down_free) {
                        std::uint64_t nc = set_plug(code, j, L);
                        nc = set_plug(nc, j + 1, 0);
                        add(nxt, nc, cnt);
                    }
                    if (right_free) {
                        std::uint64_t nc = set_plug(code, j, 0);
                        nc = set_plug(nc, j + 1, L);
                        add(nxt, nc, cnt);
                    }
                } else if (L == 0 && U != 0) {
                    // shape 2 -- extend the up path, carrying its bracket.
                    if (down_free) {
                        std::uint64_t nc = set_plug(code, j, U);
                        nc = set_plug(nc, j + 1, 0);
                        add(nxt, nc, cnt);
                    }
                    if (right_free) {
                        std::uint64_t nc = set_plug(code, j, 0);
                        nc = set_plug(nc, j + 1, U);
                        add(nxt, nc, cnt);
                    }
                } else if (L == 1 && U == 1) {
                    // shape 3 -- two opens merge. flip U's partner ")" to "(".
                    const int b = match_forward(code, j + 1);
                    std::uint64_t nc = set_plug(code, j, 0);
                    nc = set_plug(nc, j + 1, 0);
                    nc = set_plug(nc, b, 1);
                    add(nxt, nc, cnt);
                } else if (L == 2 && U == 2) {
                    // shape 3 -- two closes merge. flip L's partner "(" to ")".
                    const int b = match_backward(code, j);
                    std::uint64_t nc = set_plug(code, j, 0);
                    nc = set_plug(nc, j + 1, 0);
                    nc = set_plug(nc, b, 2);
                    add(nxt, nc, cnt);
                } else if (L == 2 && U == 1) {
                    // shape 3 -- ")" then "(", already ordered. drop both, no flip.
                    std::uint64_t nc = set_plug(code, j, 0);
                    nc = set_plug(nc, j + 1, 0);
                    add(nxt, nc, cnt);
                } else {
                    // shape 4 -- L="(", U=")" seal the same component. a CLOSE.
                    // legal only at the last free cell, and only if it empties the
                    // profile. anything else is a forbidden early sub-loop.
                    std::uint64_t nc = set_plug(code, j, 0);
                    nc = set_plug(nc, j + 1, 0);
                    if (last_cell && nc == 0) answer = (answer + cnt) % MOD;
                }
            }

            cur = std::move(nxt);
        }

        // step to the next row -- shift every profile one slot right, slot 0
        // becomes empty, slot m (forced empty past the last column) falls off.
        std::unordered_map<std::uint64_t, std::int64_t> shifted;
        shifted.reserve(cur.size() * 2);
        for (const auto& [code, cnt] : cur)
            add(shifted, (code << 2) & mask, cnt);
        cur = std::move(shifted);
    }

    return answer % MOD;
}

}  // namespace p0034
