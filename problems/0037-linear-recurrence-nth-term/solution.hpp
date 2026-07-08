#pragma once

#include <cstdint>
#include <vector>

namespace p0037 {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1000000007;  // prime -- every nonzero residue has an inverse.

// (a*b) mod p through a 128-bit intermediate. both factors sit under p < 2^30,
// so the product tops out near 1e18 -- that already fits int64. the __int128 is
// belt-and-suspenders: it also swallows the accumulation cases where we sum many
// products before reducing, and it costs a couple of instructions on arm64.
inline i64 mul(i64 a, i64 b) {
    return static_cast<i64>(static_cast<__int128>(a) * b % kMod);
}

// a^e mod p by square-and-multiply. e is small here (p-2 for inverses, or the
// bit width of whatever we raise) so this never dominates.
inline i64 pw(i64 a, i64 e) {
    a %= kMod;
    if (a < 0) a += kMod;
    i64 r = 1;
    while (e) {
        if (e & 1) r = mul(r, a);
        a = mul(a, a);
        e >>= 1;
    }
    return r;
}

// fermat inverse: p is prime, so a^(p-2) == a^{-1} (mod p). the whole reason the
// coefficient field is GF(p) and not the integers -- berlekamp-massey divides by
// the discrepancy, and division only exists because every nonzero residue has an
// inverse. this is the line that would break over a composite modulus.
inline i64 inv(i64 a) { return pw(a, kMod - 2); }

// --- berlekamp-massey over GF(p) ---
//
// hand it a_0..a_{n-1}, get back the shortest c_1..c_k with
//   a_i = c_1 a_{i-1} + c_2 a_{i-2} + … + c_k a_{i-k}   (mod p)   for all i >= k.
// this is the whole trick -- the recurrence is hidden, and BM re-discovers it
// from the terms alone in O(n*k).
//
// it walks the sequence left to right holding the current best LFSR C. at step i
// it predicts a_i from C and measures the discrepancy d = a_i - prediction. d==0
// means C still explains everything -- move on. d!=0 means C is wrong here, so it
// patches C by folding in a scaled, shifted copy of the *last* recurrence that
// failed (call it B, with its own old discrepancy b). the scale is d/b -- and
// d/b only exists because we're over a field. the patch cancels the error at i
// without disturbing the earlier fit, because B was itself shift-aligned to fail
// exactly one step later than it used to succeed.
//
// the length only grows when 2L <= i: below that the current C is already long
// enough to be patched in place; at or above it, a longer LFSR is forced, and the
// new length L' = i+1-L is provably minimal (the massey length bound). when L
// jumps we snapshot the pre-patch C as the new B -- that snapshot is the "last
// thing that worked up to here," which the next patch will lean on.
//
// C is stored as a_i + C[1] a_{i-1} + … + C[L] a_{i-L} = 0, i.e. C[0]=1 and the
// recurrence coefficients are the negatives of C[1..L]. we flip the sign on the
// way out so the caller gets c_j directly.
inline std::vector<i64> berlekamp_massey(const std::vector<i64>& s) {
    const int n = static_cast<int>(s.size());
    std::vector<i64> C(n, 0), B(n, 0), T;
    C[0] = B[0] = 1;
    int L = 0;    // length of the current LFSR.
    int mm = 0;   // steps since B was last refreshed -- the shift applied to B.
    i64 b = 1;    // discrepancy that stood when B was snapshotted.

    for (int i = 0; i < n; ++i) {
        ++mm;
        // discrepancy: a_i minus what C predicts from the k terms behind it.
        i64 d = s[i] % kMod;
        for (int j = 1; j <= L; ++j) d = (d + mul(C[j], s[i - j] % kMod)) % kMod;
        if (d == 0) continue;  // C already nails a_i -- nothing to fix.

        T = C;                       // snapshot before we mutate, in case L jumps.
        const i64 coef = mul(d, inv(b));  // d/b -- the field division BM needs.
        // C(x) -= (d/b) x^mm B(x): kill the error at i using the shifted old LFSR.
        for (int j = mm; j < n; ++j)
            C[j] = ((C[j] - mul(coef, B[j - mm])) % kMod + kMod) % kMod;

        if (2 * L > i) continue;     // C was long enough -- patched in place.
        L = i + 1 - L;               // forced growth to the minimal new length.
        B = T;                       // the pre-patch C becomes the fallback.
        b = d;
        mm = 0;
    }

    // C[0]=1, C[1..L] carry the negated recurrence. hand back c_j = -C[j].
    std::vector<i64> rec(L);
    for (int j = 1; j <= L; ++j) rec[j - 1] = ((-C[j]) % kMod + kMod) % kMod;
    return rec;
}

// --- kitamasa: a_N without touching a_1..a_{N-1} ---
//
// the recurrence is exactly the statement that the shift operator x (x·a_i =
// a_{i+1}) annihilates f(x) = x^k - c_1 x^{k-1} - … - c_k. so working modulo f,
// x^N collapses to a degree-<k polynomial:
//   x^N ≡ r_0 + r_1 x + … + r_{k-1} x^{k-1}   (mod f).
// apply both sides to a_0: x^N·a_0 = a_N on the left, and r_i x^i·a_0 = r_i a_i on
// the right. hence
//   a_N = r_0 a_0 + r_1 a_1 + … + r_{k-1} a_{k-1}.
// that's the payoff -- N never appears as a loop bound. we only ever raise x to
// the N by binary exponentiation, each step a poly multiply-then-reduce mod f.
// O(k^2 log N) total, against the O(N*k) of honest iteration. at N=1e18 that's
// the difference between ~60 squarings and a quintillion additions.
inline i64 kitamasa(const std::vector<i64>& rec, const std::vector<i64>& init,
                    u64 N) {
    const int k = static_cast<int>(rec.size());
    if (k == 0) return 0;  // empty recurrence == all-zero sequence.

    // multiply two degree-<k polynomials, then reduce mod f back to degree <k.
    // reduction rule, straight off f: x^k ≡ c_1 x^{k-1} + … + c_k, i.e.
    //   x^k ≡ sum_{j=0}^{k-1} rec[j] x^{k-1-j}.
    // for any term at position i>=k, x^i = x^{i-k}·x^k ≡ sum_j rec[j] x^{i-1-j},
    // strictly lower degrees -- so sweeping i from high to low terminates.
    auto mulmod = [&](const std::vector<i64>& A,
                      const std::vector<i64>& Bp) -> std::vector<i64> {
        // convolution first, accumulated wide. each cell sums <= k products, each
        // < 1e18; __int128 holds up to ~1.7e38, so no reduce inside the double
        // loop -- one reduction per cell after.
        std::vector<__int128> t(2 * k - 1, 0);
        for (int i = 0; i < k; ++i) {
            if (!A[i]) continue;
            const i64 ai = A[i];
            for (int j = 0; j < k; ++j) t[i + j] += static_cast<__int128>(ai) * Bp[j];
        }
        std::vector<i64> res(2 * k - 1);
        for (int i = 0; i < 2 * k - 1; ++i) res[i] = static_cast<i64>(t[i] % kMod);
        // fold the top half down through f. res[i] is final before it's read
        // because everything above i has already been cleared.
        for (int i = 2 * k - 2; i >= k; --i) {
            const i64 v = res[i];
            if (v) {
                for (int j = 0; j < k; ++j)
                    res[i - 1 - j] = (res[i - 1 - j] + mul(v, rec[j])) % kMod;
            }
            res[i] = 0;
        }
        res.resize(k);
        return res;
    };

    std::vector<i64> result(k, 0);  // the polynomial 1 == x^0.
    result[0] = 1;
    std::vector<i64> base(k, 0);    // x, reduced mod f.
    if (k == 1)
        base[0] = rec[0] % kMod;    // k==1: x ≡ c_1 already, degree 1 collapses.
    else
        base[1] = 1;                // x itself has degree 1 < k -- no reduction.

    for (u64 e = N; e; e >>= 1) {
        if (e & 1) result = mulmod(result, base);
        base = mulmod(base, base);
    }

    // a_N = sum r_i a_i. init must hold a_0..a_{k-1}.
    i64 ans = 0;
    for (int i = 0; i < k; ++i) ans = (ans + mul(result[i], init[i])) % kMod;
    return ans;
}

// end to end: recover the recurrence, then either read a_N off the prefix or leap
// to it with kitamasa.
//
// lower bound, stated honestly: you must read all m terms, so any correct solver
// is Omega(m). BM adds O(m*k) to find the order-k law; kitamasa adds
// O(k^2 log N) to reach term N. the naive alternative -- iterate a_i = sum c_j
// a_{i-j} up to N -- is O(N*k), and at N=1e18 that is 1e18 additions, dead. the
// log N in kitamasa is the entire reason this problem is tractable.
inline i64 solve(u64 N, const std::vector<i64>& a) {
    const std::vector<i64> rec = berlekamp_massey(a);
    const int k = static_cast<int>(rec.size());
    if (k == 0) return 0;                      // all-zero sequence.
    if (N < a.size()) return a[static_cast<std::size_t>(N)];  // already have it.
    const std::vector<i64> init(a.begin(), a.begin() + k);    // a_0..a_{k-1}.
    return kitamasa(rec, init, N);
}

}  // namespace p0037
