#pragma once

#include <cstdint>
#include <vector>

namespace p0037 {

using i64 = std::int64_t;
using u64 = std::uint64_t;

// same prime as the solver -- shared value, shared no logic.
inline constexpr i64 kModRef = 1000000007;

// two oracles, both dumb, both obviously correct, both handed the TRUE recurrence
// the test generated. the solution never sees that recurrence -- BM has to
// re-derive it from the terms. so agreement here is a real check: the recovered
// law and the leap-to-N have to reproduce a sequence built by an independent
// method that shares no code with either.

// oracle a -- naive_iterate. the definition of the sequence, run forward. only
// good for N a solver could actually reach in a loop (<= ~1e5 in the tests). this
// is ground truth: no cleverness to be wrong about.
//
// rec = c_1..c_k, init = a_0..a_{k-1}. a_i = sum_{j=1..k} c_j a_{i-j}.
inline i64 naive_iterate(const std::vector<i64>& rec,
                         const std::vector<i64>& init, u64 N) {
    const int k = static_cast<int>(rec.size());
    if (k == 0) return 0;  // no recurrence -- all zeros by construction.
    if (N < init.size()) return init[static_cast<std::size_t>(N)];

    std::vector<i64> seq(init);
    seq.reserve(static_cast<std::size_t>(N) + 1);
    for (u64 i = init.size(); i <= N; ++i) {
        __int128 acc = 0;  // k products under 1e18 each -- wide accumulate, one mod.
        for (int j = 0; j < k; ++j)
            acc += static_cast<__int128>(rec[j]) * seq[i - 1 - static_cast<u64>(j)];
        seq.push_back(static_cast<i64>(acc % kModRef));
    }
    return seq[static_cast<std::size_t>(N)];
}

// oracle b -- companion_matrix_power. a completely different road to a_N: no
// polynomial reduction anywhere. build the k x k companion matrix C, raise it to
// the N by fast exponentiation, read a_N off the first row. valid for ANY N,
// including 1e18, which is exactly where kitamasa's poly-mod arithmetic is most
// likely to hide a bug -- so this is the oracle that guards the hard case.
//
// C advances the state window: with v_i = (a_i, a_{i+1}, …, a_{i+k-1})^T,
//   v_{i+1} = C v_i,   C = [ 0 1 0 … 0 ]
//                          [ 0 0 1 … 0 ]
//                          [   …       ]
//                          [ 0 0 0 … 1 ]
//                          [ c_k c_{k-1} … c_1 ]
// the last row is the recurrence read backwards, so its dot with v_i is
// a_{i+k}. therefore v_N = C^N v_0 and a_N = (C^N v_0)_0 = row0(C^N) · init.
class CompanionOracle {
public:
    static constexpr i64 kModRef = 1000000007;

    // returns a_N for the recurrence rec with seeds init, any N up to 1e18.
    static i64 nth(const std::vector<i64>& rec, const std::vector<i64>& init,
                   u64 N) {
        const int k = static_cast<int>(rec.size());
        if (k == 0) return 0;

        Mat C(k, std::vector<i64>(k, 0));
        for (int i = 0; i + 1 < k; ++i) C[i][i + 1] = 1;  // shift block.
        for (int j = 0; j < k; ++j) C[k - 1][j] = rec[k - 1 - j];  // c_k..c_1.

        const Mat P = matpow(C, N, k);
        i64 ans = 0;  // a_N = row 0 of C^N, dotted with the seeds.
        for (int j = 0; j < k; ++j)
            ans = (ans + mul(P[0][j], init[j])) % kModRef;
        return ans;
    }

private:
    using Mat = std::vector<std::vector<i64>>;

    static i64 mul(i64 a, i64 b) {
        return static_cast<i64>(static_cast<__int128>(a) * b % kModRef);
    }

    static Mat matmul(const Mat& A, const Mat& B, int k) {
        Mat R(k, std::vector<i64>(k, 0));
        for (int i = 0; i < k; ++i)
            for (int t = 0; t < k; ++t) {
                if (!A[i][t]) continue;
                const i64 a = A[i][t];
                for (int j = 0; j < k; ++j)
                    R[i][j] = (R[i][j] + mul(a, B[t][j])) % kModRef;
            }
        return R;
    }

    static Mat matpow(Mat base, u64 e, int k) {
        Mat r(k, std::vector<i64>(k, 0));
        for (int i = 0; i < k; ++i) r[i][i] = 1;  // identity: C^0.
        while (e) {
            if (e & 1) r = matmul(r, base, k);
            base = matmul(base, base, k);
            e >>= 1;
        }
        return r;
    }
};

// thin wrapper so tests read the same on both oracles.
inline i64 companion_matrix_power(const std::vector<i64>& rec,
                                  const std::vector<i64>& init, u64 N) {
    return CompanionOracle::nth(rec, init, N);
}

}  // namespace p0037
