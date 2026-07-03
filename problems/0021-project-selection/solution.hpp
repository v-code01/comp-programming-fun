#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace p0021 {

// project selection -- the maximum-weight closure problem, solved by min-cut.
//
// systems framing: you're deciding which build targets / features / projects to
// ship. each one has a net profit p_i (a target that pulls its weight is +, one
// that only costs is -). selecting a target drags in mandatory prerequisites:
// relation i->j means "ship i and you must ship j too." pick the subset that
// maximizes net profit and is closed under those prerequisites. ship nothing and
// you net zero, so the answer never dips below 0.
//
// this is textbook. a subset S is a CLOSURE when no prerequisite edge leaves it:
// i in S and i->j forces j in S. maximum-weight closure reduces to min-cut, and
// min-cut is max-flow by the max-flow/min-cut theorem (Ford & Fulkerson, 1956).
//
// ---- the reduction, and why it's exact ----
//
// build a flow network on n+2 nodes: source s, sink t, one node per project.
//   p_i > 0:  s -> i  with capacity  p_i
//   p_i < 0:  i -> t  with capacity  |p_i|
//   i -> j (prereq): i -> j with capacity INF   -- never worth cutting.
//
// take any finite s-t cut and let A be the source side. an INF edge crossing the
// cut would blow the capacity past anything finite, so a finite cut never cuts
// one -- meaning if i in A and i->j exists, j is in A too. so A minus s is a
// valid closure. the cut only severs s->i edges (for positive i left out of A)
// and i->t edges (for negative i pulled into A). its capacity is
//   cut(A) = sum_{i not in A, p_i>0} p_i  +  sum_{i in A, p_i<0} |p_i|.
//
// let P = sum of all positive p_i. the profit of the closure S = A minus s is
//   profit(S) = sum_{i in S} p_i
//             = sum_{i in S, p_i>0} p_i  -  sum_{i in S, p_i<0} |p_i|
//             = (P - sum_{i not in S, p_i>0} p_i) - sum_{i in S, p_i<0} |p_i|
//             = P - cut(A).
// every finite cut is a closure and every closure is a finite cut, term for term.
// so maximizing profit is minimizing the cut -- and min-cut = max-flow. therefore
//   max profit = P - maxflow(s, t).
// the empty closure is the cut A = {s}: it severs every positive s-edge, capacity
// P, profit 0. so maxflow <= P and the answer is >= 0, as required.
//
// INF must exceed any finite cut, i.e. exceed P. worst case P = 5000 * 1e9 = 5e12,
// so 1e18 clears it by six orders of magnitude and still leaves int64 headroom:
// flow never saturates an INF edge, so no capacity ever grows past P + a bit.
inline constexpr std::int64_t kInf = 1'000'000'000'000'000'000LL;  // 1e18 > 5e12.

// Dinic's algorithm. O(V^2 * E) worst case, far faster on these graphs in
// practice -- BFS carves the residual into levels, DFS pushes a blocking flow
// along shortest augmenting paths, repeat until t is unreachable. capacities are
// int64 because both the profits and the running flow live at the 1e12..1e18
// scale.
class Dinic {
public:
    explicit Dinic(int num_nodes)
        : head_(static_cast<std::size_t>(num_nodes), -1),
          level_(static_cast<std::size_t>(num_nodes)),
          iter_(static_cast<std::size_t>(num_nodes)),
          n_(num_nodes) {}

    // a directed edge with capacity cap, plus its residual twin (capacity 0).
    // storing both in one flat array means edge e's twin is e^1 -- push flow on
    // one, pull it back on the other, no pointer chasing.
    void add_edge(int from, int to, std::int64_t cap) {
        to_.push_back(to);
        cap_.push_back(cap);
        next_.push_back(head_[static_cast<std::size_t>(from)]);
        head_[static_cast<std::size_t>(from)] = static_cast<int>(to_.size()) - 1;

        to_.push_back(from);
        cap_.push_back(0);
        next_.push_back(head_[static_cast<std::size_t>(to)]);
        head_[static_cast<std::size_t>(to)] = static_cast<int>(to_.size()) - 1;
    }

    std::int64_t max_flow(int s, int t) {
        std::int64_t flow = 0;
        while (bfs(s, t)) {
            // reset the per-node edge cursor -- each phase walks every adjacency
            // list once, that's what bounds the DFS work.
            for (int v = 0; v < n_; ++v) {
                iter_[static_cast<std::size_t>(v)] = head_[static_cast<std::size_t>(v)];
            }
            std::int64_t f;
            while ((f = dfs(s, t, kInf)) > 0) flow += f;
        }
        return flow;
    }

private:
    // level graph: BFS distance from s in the residual network. edges only go
    // level d -> level d+1 during the DFS, which is what makes augmenting paths
    // shortest and bounds the phase count to O(V).
    bool bfs(int s, int t) {
        for (auto& l : level_) l = -1;
        level_[static_cast<std::size_t>(s)] = 0;
        // ring-free BFS on a reused vector -- no per-call allocation in the hot
        // loop, we just index a growing frontier.
        bfs_q_.clear();
        bfs_q_.push_back(s);
        for (std::size_t qi = 0; qi < bfs_q_.size(); ++qi) {
            const int v = bfs_q_[qi];
            for (int e = head_[static_cast<std::size_t>(v)]; e != -1;
                 e = next_[static_cast<std::size_t>(e)]) {
                const int u = to_[static_cast<std::size_t>(e)];
                if (cap_[static_cast<std::size_t>(e)] > 0 &&
                    level_[static_cast<std::size_t>(u)] < 0) {
                    level_[static_cast<std::size_t>(u)] =
                        level_[static_cast<std::size_t>(v)] + 1;
                    bfs_q_.push_back(u);
                }
            }
        }
        return level_[static_cast<std::size_t>(t)] >= 0;
    }

    // push up to f units from v to t along level-increasing residual edges.
    // iter_ remembers where each node left off, so a dead edge is skipped once
    // per phase, never revisited. recursion depth is bounded by the number of
    // levels (<= V), and each level advances by one, so the stack stays under V
    // frames -- safe for n up to 5000 on any real stack.
    std::int64_t dfs(int v, int t, std::int64_t f) {
        if (v == t) return f;
        for (int& e = iter_[static_cast<std::size_t>(v)]; e != -1;
             e = next_[static_cast<std::size_t>(e)]) {
            const int u = to_[static_cast<std::size_t>(e)];
            if (cap_[static_cast<std::size_t>(e)] > 0 &&
                level_[static_cast<std::size_t>(u)] ==
                    level_[static_cast<std::size_t>(v)] + 1) {
                const std::int64_t d =
                    dfs(u, t, f < cap_[static_cast<std::size_t>(e)]
                                  ? f
                                  : cap_[static_cast<std::size_t>(e)]);
                if (d > 0) {
                    cap_[static_cast<std::size_t>(e)] -= d;       // spend forward.
                    cap_[static_cast<std::size_t>(e ^ 1)] += d;   // credit the twin.
                    return d;
                }
            }
        }
        return 0;
    }

    // flat forward-star adjacency. to_/cap_/next_ index edges; twin of e is e^1.
    std::vector<int> head_;
    std::vector<int> to_;
    std::vector<std::int64_t> cap_;
    std::vector<int> next_;
    std::vector<int> level_;
    std::vector<int> iter_;
    std::vector<int> bfs_q_;
    int n_;
};

// the whole solve. profits are 0-indexed (project k is profits[k]); prereqs are
// 0-indexed pairs (i, j) meaning "selecting i requires j." returns the max net
// profit of any prerequisite-closed subset, always >= 0.
inline std::int64_t max_profit(const std::vector<std::int64_t>& profits,
                               const std::vector<std::pair<int, int>>& prereqs) {
    const int n = static_cast<int>(profits.size());
    // node layout: 0 = source, 1..n = projects, n+1 = sink.
    const int s = 0;
    const int t = n + 1;
    Dinic din(n + 2);

    std::int64_t positive_sum = 0;  // this is P -- the sum of all positive profits.
    for (int i = 0; i < n; ++i) {
        if (profits[static_cast<std::size_t>(i)] > 0) {
            positive_sum += profits[static_cast<std::size_t>(i)];
            din.add_edge(s, i + 1, profits[static_cast<std::size_t>(i)]);
        } else if (profits[static_cast<std::size_t>(i)] < 0) {
            din.add_edge(i + 1, t, -profits[static_cast<std::size_t>(i)]);
        }
        // p_i == 0 attaches to neither side -- it's free to include or drop and
        // never changes the profit or the cut.
    }
    // prerequisite edges get INF capacity: the min cut never severs one, so the
    // source side is forced to stay closed.
    for (const auto& pr : prereqs) {
        din.add_edge(pr.first + 1, pr.second + 1, kInf);
    }

    // max profit = P - min-cut = P - max-flow.
    return positive_sum - din.max_flow(s, t);
}

}  // namespace p0021
