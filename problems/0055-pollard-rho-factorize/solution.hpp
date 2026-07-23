#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace p0055 {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// COMPLEXITY. reading is Omega(q). one number costs ~O(n^{1/4} log n): pollard's
// rho finds a factor in expected O(n^{1/4}) iterations (birthday bound on the
// residues of the smallest prime factor), each iteration a handful of mulmods,
// and miller-rabin is O(log n) modpows to certify primality. trial division is
// O(sqrt n) = 1e9 at n = 1e18 -- dead on arrival. this is the reason rho exists.

// multiply two residues mod m, both < m <= 1e18. the product tops out near 1e36,
// and an unsigned __int128 reaches 3.4e38 -- so it never wraps. that's the whole
// trick behind 64-bit factoring: one 128-bit multiply, one 128-bit remainder,
// no overflow. montgomery multiplication would drop the divide for a real
// speedup, but a single 128-bit % is correct and simplest -- reach for
// montgomery only when the profile asks for it.
inline u64 mulmod(u64 a, u64 b, u64 m) {
    return static_cast<u64>(static_cast<u128>(a) * b % m);
}

// a^e mod m by square-and-multiply. O(log e) mulmods.
inline u64 modpow(u64 a, u64 e, u64 m) {
    u64 r = 1 % m;  // m == 1 folds to 0; witness paths guard m > 37 anyway.
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

// deterministic miller-rabin over the whole 64-bit range. write n-1 = d * 2^s.
// a base a is a liar for a prime only if a^d == 1 or a^(d*2^r) == -1 for some
// r < s. the twelve bases {2,3,5,7,11,13,17,19,23,29,31,37} -- every prime below
// 41 -- have NO composite that survives all of them below 3.317e24 (jaeschke,
// tightened by sorenson-webster). 1e18 sits far under that bound, so this is
// EXACT here, not probabilistic: no false positive, no false negative. small n
// and even n fall out of the trial loop up front.
inline bool is_prime(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull,
                  29ull, 31ull, 37ull}) {
        if (n % p == 0) return n == p;  // divisible by a small prime -> that's it.
    }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }
    for (u64 a : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull,
                  29ull, 31ull, 37ull}) {
        u64 x = modpow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = mulmod(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;  // a witnessed compositeness -- done.
    }
    return true;
}

inline u64 gcd(u64 a, u64 b) {
    while (b) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;  // gcd(0, n) == n -- the batched-collision signal below leans on it.
}

// the rho map. f(x) = x^2 + c mod n, a cheap pseudo-random walk on Z/n.
inline u64 rho_f(u64 x, u64 c, u64 n) { return (mulmod(x, x, n) + c) % n; }

// pollard's rho, brent's variant, gcd batched. iterating f mod n eventually
// repeats a residue mod some prime factor p, and by the birthday bound that
// collision lands after ~sqrt(p) <= n^{1/4} steps -- ~3e4 for a 1e18 semiprime,
// not the 1e9 that trial division would pay. brent swaps floyd's two-pointer
// chase for power-of-two teleports (the r doubling), and we fold the step
// differences into one running product q so it's a single gcd per ~128 steps
// instead of one per step. the map has fixed points that stall or collide the
// whole batch at once -- when that shows up as gcd == n, backtrack step by step,
// and if the backtrack is still degenerate bump c and rewalk. only ever called
// on composite n, so a nontrivial divisor exists and some c finds it.
inline u64 pollard_rho(u64 n) {
    if ((n & 1) == 0) return 2;
    const u64 batch = 128;
    for (u64 c = 1;; ++c) {
        u64 x = 2, y = 2, ys = 2, d = 1, q = 1, r = 1;
        do {
            x = y;
            for (u64 i = 0; i < r; ++i) y = rho_f(y, c, n);
            for (u64 k = 0; k < r && d == 1; k += batch) {
                ys = y;
                u64 lim = std::min<u64>(batch, r - k);
                for (u64 i = 0; i < lim; ++i) {
                    y = rho_f(y, c, n);
                    u64 diff = x > y ? x - y : y - x;
                    q = mulmod(q, diff, n);  // zero diff zeroes q -> gcd hits n.
                }
                d = gcd(q, n);
            }
            r <<= 1;
        } while (d == 1);
        if (d != n) return d;
        // the batch swallowed the factor -- replay the last run one step at a
        // time to isolate it.
        do {
            ys = rho_f(ys, c, n);
            u64 diff = x > ys ? x - ys : ys - x;
            d = gcd(diff, n);
        } while (d == 1);
        if (d != n) return d;
        // c was unlucky (cycle closed with no proper factor). bump it, rewalk.
    }
}

// factor n into primes. prime by miller-rabin -> emit it; else split off a
// nontrivial divisor with rho and recurse on both halves. n == 1 contributes
// nothing.
inline void factor_rec(u64 n, std::vector<u64>& out) {
    if (n == 1) return;
    if (is_prime(n)) {
        out.push_back(n);
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

// the prime factors of n in non-decreasing order, with multiplicity. n == 1
// gives the empty list.
inline std::vector<u64> factorize(u64 n) {
    std::vector<u64> out;
    factor_rec(n, out);
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace p0055
