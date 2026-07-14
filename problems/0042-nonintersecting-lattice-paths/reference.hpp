#pragma once

#include <cstdint>
#include <vector>

#include "solution.hpp"  // for Pt only. no binomials, no determinant, no field.

namespace p0042 {

// the dumb oracle. it shares NOTHING with LGV -- no binomial, no matrix, no
// cancellation argument, no permutation. it just builds the paths and looks.
//
// for each i: walk every monotone path from a_i to b_i with a DFS over {right,
// up}. that's the complete list, by construction -- at every lattice point the
// only two moves are right and up, and we only take one when it doesn't overshoot
// the destination.
//
// then take the cartesian product of the k lists and count the tuples whose paths
// are pairwise vertex-disjoint. disjointness is a set question, so a path is
// carried as its set of lattice points, packed one bit per point.
//
// exponential, obviously. C(dx+dy, dx) paths per pair, then a k-fold product. it
// is the ground truth for tiny inputs and nothing else. spans of ~4 and k <= 3.
class ReferenceSolver {
public:
    // the grid the oracle will look at: 0..7 in both axes, so a point is one bit
    // of a u64 at index y*8+x. anything outside and the packing would alias two
    // different points onto one bit -- so refuse loudly with -1 instead of
    // quietly returning a wrong count.
    static constexpr int kMaxCoord = 7;

    // families with path i running a_i -> b_perm[i], all k paths pairwise
    // vertex-disjoint. pass the identity for the problem as stated; pass anything
    // else to interrogate the compatibility condition itself.
    std::int64_t solve(const std::vector<Pt>& a, const std::vector<Pt>& b,
                       const std::vector<int>& perm) {
        const int k = static_cast<int>(a.size());
        for (const Pt& p : a)
            if (!in_grid(p)) return -1;
        for (const Pt& p : b)
            if (!in_grid(p)) return -1;

        lists_.assign(static_cast<std::size_t>(k), {});
        for (int i = 0; i < k; ++i) {
            const Pt s = a[static_cast<std::size_t>(i)];
            const Pt t = b[static_cast<std::size_t>(perm[static_cast<std::size_t>(i)])];
            if (t.x < s.x || t.y < s.y) continue;  // unreachable -- empty list, zero families.
            walk(s, t, 0, lists_[static_cast<std::size_t>(i)]);
        }

        count_ = 0;
        pick(0, k, 0);
        return count_;
    }

    // the problem's pairing: path i goes a_i -> b_i.
    std::int64_t solve(const std::vector<Pt>& a, const std::vector<Pt>& b) {
        std::vector<int> id(a.size());
        for (std::size_t i = 0; i < id.size(); ++i) id[i] = static_cast<int>(i);
        return solve(a, b, id);
    }

private:
    static bool in_grid(Pt p) {
        return p.x >= 0 && p.y >= 0 && p.x <= kMaxCoord && p.y <= kMaxCoord;
    }
    static std::uint64_t bit(Pt p) {
        return 1ULL << (static_cast<unsigned>(p.y) * 8u + static_cast<unsigned>(p.x));
    }

    // every monotone path from cur to dst, as a point set. cur is always on the
    // path, hence the |= before anything else.
    void walk(Pt cur, Pt dst, std::uint64_t mask, std::vector<std::uint64_t>& out) {
        mask |= bit(cur);
        if (cur == dst) {
            out.push_back(mask);
            return;
        }
        if (cur.x < dst.x) walk(Pt{cur.x + 1, cur.y}, dst, mask, out);  // right.
        if (cur.y < dst.y) walk(Pt{cur.x, cur.y + 1}, dst, mask, out);  // up.
    }

    // cartesian product, with the one prune that keeps it finishable: a partial
    // tuple that already overlaps can never become disjoint, so drop it now.
    // pruning changes nothing about what's counted -- only how long it takes.
    void pick(int i, int k, std::uint64_t used) {
        if (i == k) {
            ++count_;
            return;
        }
        for (std::uint64_t m : lists_[static_cast<std::size_t>(i)]) {
            if (m & used) continue;  // shares a lattice point with an earlier path.
            pick(i + 1, k, used | m);
        }
    }

    std::vector<std::vector<std::uint64_t>> lists_;
    std::int64_t count_ = 0;
};

inline std::int64_t solve_ref(const std::vector<Pt>& a, const std::vector<Pt>& b) {
    ReferenceSolver s;
    return s.solve(a, b);
}

inline std::int64_t solve_ref(const std::vector<Pt>& a, const std::vector<Pt>& b,
                              const std::vector<int>& perm) {
    ReferenceSolver s;
    return s.solve(a, b, perm);
}

}  // namespace p0042
