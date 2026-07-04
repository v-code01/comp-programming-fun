#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace p0029 {

// text T and n patterns P_1..P_n (repeats allowed). count, over every pattern,
// how many times it sits inside T as a substring -- overlaps counted, a pattern
// listed twice counted twice -- and sum it. (|T|, n, sum|P| all <= 1e5, lowercase.
// the answer can reach ~|T| * n, so it's an int64. a text of 1e5 a's against 1e5
// copies of "a" is 1e10 occurrences -- an int overflows, this doesn't.)
//
// the naive answer scans T once per pattern: O(|T| * sum|P|), quadratic in the
// worst case. aho-corasick reads T exactly once no matter how many patterns
// there are. that's the whole trade.
//
// LOWER BOUND. you have to read T and every pattern -- Omega(|T| + sum|P|), each
// character looked at once. the build touches every trie edge across the whole
// alphabet, O(sum|P| * 26); the text run is a single pass, O(|T|). so the floor
// and the ceiling meet up to the alphabet constant -- this is optimal.
//
// THE AUTOMATON. drop all patterns into a trie. cnt_end_[node] = how many
// patterns end exactly here (a repeat bumps it again -- that's how a listed-twice
// pattern gets counted twice). then two things make the text run O(1) per char:
//
//   FAIL LINKS. fail_[v] points at the longest proper suffix of v's string that
//   is also a trie node. build them by BFS -- when you reach v via edge c from
//   its parent p, the failure of v is "follow p's failure, then take edge c",
//   i.e. go_[fail_[p]][c]. p is shallower than v, so its failure and its goto
//   are already final by the time v is dequeued. the root and its children fail
//   to the root.
//
//   GOTO. instead of walking the fail chain at query time, precompute the full
//   transition go_[v][c] for every state and letter. a present trie edge stays
//   put; a missing one inherits go_[fail_[v]][c]. after BFS, go_ is a complete
//   DFA -- one array read advances the text, no chain-walk, no branch.
//
// THE SUFFIX-LINK COUNT SUM. when the text sits in state v, every pattern that
// ends at v OR at any node up its fail chain has just occurred ending here. that
// chain sum is fixed per state, so precompute it once:
//
//   sum_end_[v] = cnt_end_[v] + sum_end_[fail_[v]]
//
// fail_[v] is strictly shallower, so in BFS order it's already done -- one pass
// fills the table. now feeding T is: cur = go_[cur][c]; total += sum_end_[cur].
// no fail-chain walk at query time -- that walk is where the naive suffix-link
// version turns O(|T| * max_depth), and this is what kills it.
class AhoCorasick {
public:
    static constexpr int kAlpha = 26;  // lowercase 'a'..'z'.

    // build the DFA over the patterns. clears any prior build -- one instance,
    // reused across tests, must not carry state between builds.
    void build(const std::vector<std::string>& patterns) {
        reset();
        for (const std::string& p : patterns) insert(p);
        build_links();
    }

    // feed T through the DFA and sum the per-state chain counts. this is the
    // measured hot loop: one goto read and one add per character.
    std::int64_t count(const std::string& text) const {
        std::int64_t total = 0;
        int cur = 0;  // start at the root.
        for (unsigned char ch : text) {
            const int c = static_cast<int>(ch) - 'a';
            cur = go_[static_cast<std::size_t>(cur) * kAlpha +
                      static_cast<std::size_t>(c)];
            total += sum_end_[static_cast<std::size_t>(cur)];
        }
        return total;
    }

    int states() const { return node_count_; }

private:
    // go_ is a flat states*26 table: go_[v*26 + c] is the DFA transition. flat
    // over a vector<array> keeps the whole automaton in one contiguous block --
    // the text run streams it, so cache locality on this array is the run time.
    std::vector<int> go_;
    std::vector<int> fail_;
    std::vector<int> cnt_end_;         // patterns ending exactly at a node.
    std::vector<std::int64_t> sum_end_;  // that count summed up the fail chain.
    int node_count_ = 0;

    void reset() {
        go_.clear();
        fail_.clear();
        cnt_end_.clear();
        sum_end_.clear();
        node_count_ = 0;
        add_node();  // node 0 is the root.
    }

    // append one fresh node with no edges. -1 marks "no trie child" -- it stays
    // -1 through insertion and gets overwritten with a real goto during BFS. 0
    // can't stand in for "empty" here because 0 is the root itself.
    int add_node() {
        go_.insert(go_.end(), kAlpha, -1);
        fail_.push_back(0);
        cnt_end_.push_back(0);
        sum_end_.push_back(0);
        return node_count_++;
    }

    // walk the pattern from the root, minting nodes for missing edges, and mark
    // the terminal. ++ not = so a duplicate pattern raises the count again.
    void insert(const std::string& p) {
        int cur = 0;
        for (unsigned char ch : p) {
            const int c = static_cast<int>(ch) - 'a';
            const std::size_t slot =
                static_cast<std::size_t>(cur) * kAlpha + static_cast<std::size_t>(c);
            if (go_[slot] == -1) go_[slot] = add_node();
            cur = go_[slot];
        }
        ++cnt_end_[static_cast<std::size_t>(cur)];
    }

    // BFS from the root: fix fail links, complete the goto table, then roll the
    // chain-count sum down in the same visitation order.
    void build_links() {
        std::vector<int> queue;
        queue.reserve(static_cast<std::size_t>(node_count_));

        // root's children fail to the root; every absent root edge is a self-loop
        // so a mismatch at the root keeps the text sitting at the root.
        for (int c = 0; c < kAlpha; ++c) {
            const std::size_t slot = static_cast<std::size_t>(c);  // root = node 0.
            int& u = go_[slot];
            if (u == -1) {
                u = 0;
            } else {
                fail_[static_cast<std::size_t>(u)] = 0;
                queue.push_back(u);
            }
        }

        // process in BFS order. head walks the queue; children get appended, so
        // a node is always dequeued after every strictly-shallower node -- which
        // is exactly what makes fail_[u] and sum_end_ available when needed.
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const int v = queue[head];
            const std::size_t vbase = static_cast<std::size_t>(v) * kAlpha;
            const int fv = fail_[static_cast<std::size_t>(v)];
            const std::size_t fbase = static_cast<std::size_t>(fv) * kAlpha;
            for (int c = 0; c < kAlpha; ++c) {
                int& u = go_[vbase + static_cast<std::size_t>(c)];
                const int through_fail = go_[fbase + static_cast<std::size_t>(c)];
                if (u == -1) {
                    // no trie edge here -- inherit the failure state's transition.
                    u = through_fail;
                } else {
                    // real edge: its failure is where the parent's failure goes on c.
                    fail_[static_cast<std::size_t>(u)] = through_fail;
                    queue.push_back(u);
                }
            }
        }

        // sum_end_[v] = cnt_end_[v] + sum_end_[fail_[v]]. root first (its fail is
        // itself and it holds no terminal), then the queue in BFS order.
        sum_end_[0] = cnt_end_[0];
        for (std::size_t i = 0; i < queue.size(); ++i) {
            const int v = queue[i];
            sum_end_[static_cast<std::size_t>(v)] =
                cnt_end_[static_cast<std::size_t>(v)] +
                sum_end_[static_cast<std::size_t>(fail_[static_cast<std::size_t>(v)])];
        }
    }
};

// total occurrences of all patterns in text. builds the automaton, runs the text
// once. O(sum|P| * 26 + |T|) time, O(states * 26) space.
inline std::int64_t total_occurrences(const std::string& text,
                                      const std::vector<std::string>& patterns) {
    AhoCorasick ac;
    ac.build(patterns);
    return ac.count(text);
}

}  // namespace p0029
