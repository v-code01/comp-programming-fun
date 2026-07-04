#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace p0031 {

// wildcard string matching. text T is lowercase, pattern P is lowercase with
// '?' meaning "any single char". find every 0-indexed start i where P lines up
// on T[i .. i+m-1], '?' matching anything.
//
// the naive check is O(n*m) -- slide P over T, compare m chars per stop. this
// does O((n+m) log(n+m)) with two FFTs. here's the whole trick.
//
// give each letter a value: a..z -> 1..26. text value T[t], pattern value P[j].
// weight w[j] = 0 when P[j]=='?', else 1. then position i matches iff
//
//     M(i) = sum_j  w[j] * (T[i+j] - P[j])^2  == 0.
//
// every term is w>=0 times a square >=0, so the sum is zero only when EVERY term
// is zero -- meaning for each j either P[j] is a wildcard (w=0) or T[i+j]==P[j].
// that's the match condition, exactly. no approximation in the model, only in
// the arithmetic (see PRECISION below).
//
// expand the square and split by which index the sum runs over:
//
//     M(i) = sum_j w[j] T[i+j]^2  -  2 sum_j (w[j] P[j]) T[i+j]  +  sum_j w[j] P[j]^2.
//
// the third piece runs over pattern indices j only -- it does NOT depend on i.
// P doesn't slide. so it's a single scalar
//
//     Csum = sum_j w[j] P[j]^2,
//
// computed once, not a convolution. the textbook framing calls this "three
// correlations"; the P^2 one degenerates to a constant because the pattern is
// fixed. two correlations actually vary with i:
//
//     C1(i) = sum_j A[j] * T2[i+j]      A[j] = w[j],           T2[t] = T[t]^2
//     C2(i) = sum_j B[j] * T[i+j]       B[j] = w[j] P[j].
//
// a correlation sum_j A[j] X[i+j] becomes a convolution once you reverse A:
// let Arev[k] = A[m-1-k], then conv(Arev, X)[i+m-1] = sum_j A[j] X[i+j]. so
//
//     M(i) = conv(Arev, T2)[i+m-1] - 2*conv(Brev, T)[i+m-1] + Csum,
//
// two FFT convolutions and a scalar. match iff M(i) rounds to 0.
//
// LOWER BOUND, honestly. you must read T and P -- Omega(n+m). the two FFTs are
// O((n+m) log(n+m)). that log factor is the whole price over the Omega(n+m)
// floor, and it buys you the beating of the naive O(n*m): at n=m=2e5 that's
// 4e10 char compares vs ~4e5*19 ~ 8e6 butterflies. not close.
//
// FFT vs NTT. picked FFT (double). the values are small and bounded, so the
// rounding is provably safe with room to spare (below), and doubles skip the
// modular-inverse and CRT machinery an exact NTT needs. NTT is the exact
// alternative -- if the bounds grew until eps*V crept toward 0.5 you'd switch.
// they don't, so we don't.
//
// PRECISION -- the one place this problem bites. M(i) is a sum of at most m
// nonnegative terms, each w*(T-P)^2 <= 1 * 26^2 = 676, so 0 <= M(i) <= 676*m.
// at m=2e5 that caps M at ~1.35e8, and it's an integer. the FFT computes it in
// doubles, so the question is only: is the accumulated float error < 0.5, so
// that llround recovers the exact integer?
//   the convolutions carry the magnitude. conv(Arev,T2) tops out where the
// transforms multiply: |FFT(Arev)| <= sum|A| <= m ~ 2e5, |FFT(T2)| <= sum(T2)
// <= 676n ~ 1.35e8, so the pointwise product sits at ~2.7e13 -- comfortably
// inside a double's 2^53 ~ 9e15 exact-integer reach, no mantissa loss in the
// hot step. (this is WHY the two real inputs get unpacked into separate
// transforms before multiplying, rather than the square-in-place packing --
// squaring would push a 1.35e8 DC term to 1.8e16 and start shaving bits.)
//   the standard FFT error model gives per-output absolute error on the order
// of  eps * log2(N) * Vmax, with eps=2.2e-16, N<=2^19 so log2 N ~ 19, and
// Vmax ~ 3e8 (the 2*C2 term doubles one contribution). that's ~1.3e-6. five-ish
// orders of magnitude under 0.5. round to nearest, compare == 0, done. the
// large-magnitude tests (all-'z' text and pattern, full 2e5 length) exercise
// this exact worst case against the direct oracle.

// ---- a small complex, hand-rolled ---------------------------------------
// std::complex<double>'s operator* guards against inf/nan on every multiply --
// dead weight in an FFT inner loop. this doesn't.
struct Cd {
    double re, im;
};

inline Cd operator+(Cd a, Cd b) { return {a.re + b.re, a.im + b.im}; }
inline Cd operator-(Cd a, Cd b) { return {a.re - b.re, a.im - b.im}; }
inline Cd operator*(Cd a, Cd b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

// iterative in-place radix-2 forward DFT. n must be a power of two.
// the roots table is built once and grown across calls (static): each root is a
// product of at most log2(n) exact factors, so its error stays ~log2(n)*eps
// instead of the len/2*eps you'd accumulate multiplying a running twiddle down
// a segment -- that difference is the margin the precision argument leans on.
// inverse transform is done by the caller via conjugation, so this stays
// forward-only and single-purpose.
inline void fft(std::vector<Cd>& a) {
    const int n = static_cast<int>(a.size());
    if (n == 1) return;  // size-1 DFT is the identity; the shifts below would UB.

    static std::vector<Cd> rt = {{1, 0}, {1, 0}};
    static int built = 2;  // rt is valid up to this size; extend as n grows.
    if (built < n) {
        rt.resize(static_cast<std::size_t>(n));
        for (; built < n; built <<= 1) {
            const double ang = std::acos(-1.0) / built;  // pi / built.
            const Cd x{std::cos(ang), std::sin(ang)};
            for (int i = built; i < 2 * built; ++i)
                rt[static_cast<std::size_t>(i)] =
                    (i & 1) ? rt[static_cast<std::size_t>(i >> 1)] * x
                            : rt[static_cast<std::size_t>(i >> 1)];
        }
    }

    // bit-reversal permutation, derived in one pass from the half index.
    std::vector<int> rev(static_cast<std::size_t>(n));
    const int L = __builtin_ctz(static_cast<unsigned>(n));
    for (int i = 0; i < n; ++i)
        rev[static_cast<std::size_t>(i)] =
            (rev[static_cast<std::size_t>(i >> 1)] >> 1) | ((i & 1) << (L - 1));
    for (int i = 0; i < n; ++i)
        if (i < rev[static_cast<std::size_t>(i)])
            std::swap(a[static_cast<std::size_t>(i)],
                      a[static_cast<std::size_t>(rev[static_cast<std::size_t>(i)])]);

    for (int k = 1; k < n; k <<= 1)
        for (int i = 0; i < n; i += 2 * k)
            for (int j = 0; j < k; ++j) {
                const Cd z = rt[static_cast<std::size_t>(j + k)] *
                             a[static_cast<std::size_t>(i + j + k)];
                a[static_cast<std::size_t>(i + j + k)] =
                    a[static_cast<std::size_t>(i + j)] - z;
                a[static_cast<std::size_t>(i + j)] =
                    a[static_cast<std::size_t>(i + j)] + z;
            }
}

// linear convolution of two real sequences. res[k] = sum_i a[i] b[k-i].
// packs a into the real part and b into the imaginary part of one complex
// array, one forward FFT, then unpacks the two real transforms via conjugate
// symmetry -- so the pointwise product is FFT(a)*FFT(b) at ~2.7e13, not a
// squared 1.8e16. inverse via conjugation reuses the same forward routine.
inline std::vector<double> convolve(const std::vector<double>& a,
                                    const std::vector<double>& b) {
    if (a.empty() || b.empty()) return {};
    const std::size_t need = a.size() + b.size() - 1;
    std::size_t n = 1;
    while (n < need) n <<= 1;

    std::vector<Cd> in(n, Cd{0.0, 0.0});
    for (std::size_t i = 0; i < a.size(); ++i) in[i].re = a[i];
    for (std::size_t i = 0; i < b.size(); ++i) in[i].im = b[i];
    fft(in);

    // unpack. with F = FFT(a + i b), conjugate symmetry gives
    //   FA[k] = (F[k] + conj(F[n-k])) / 2
    //   FB[k] = (F[k] - conj(F[n-k])) / (2i)
    // and the transform of the convolution is FA[k]*FB[k].
    std::vector<Cd> out(n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t j = (n - k) & (n - 1);
        const Cd conj_nk{in[j].re, -in[j].im};
        const Cd fa{0.5 * (in[k].re + conj_nk.re), 0.5 * (in[k].im + conj_nk.im)};
        const Cd d = in[k] - conj_nk;      // F[k] - conj(F[n-k])
        const Cd fb{0.5 * d.im, -0.5 * d.re};  // d / (2i)
        out[k] = fa * fb;
    }

    // inverse FFT via  IFFT(x) = conj(FFT(conj(x))) / n. we only want the real
    // part, and conjugation leaves the real part untouched, so: conjugate,
    // forward-transform, divide by n, read .re.
    for (std::size_t k = 0; k < n; ++k) out[k].im = -out[k].im;
    fft(out);
    const double invn = 1.0 / static_cast<double>(n);
    std::vector<double> res(need);
    for (std::size_t i = 0; i < need; ++i) res[i] = out[i].re * invn;
    return res;
}

// every 0-indexed start i in `text` where `pat` matches, ascending. '?' in pat
// matches any char; text is plain lowercase. empty when nothing matches or when
// |pat| > |text|.
inline std::vector<int> matches(const std::string& text, const std::string& pat) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pat.size());
    std::vector<int> res;
    if (m == 0 || m > n) return res;  // no window to place the pattern in.

    // text values 1..26, and its square. lowercase only, so no wildcard here.
    std::vector<double> T(static_cast<std::size_t>(n));
    std::vector<double> T2(static_cast<std::size_t>(n));
    for (int t = 0; t < n; ++t) {
        const double v = static_cast<double>(text[static_cast<std::size_t>(t)] - 'a' + 1);
        T[static_cast<std::size_t>(t)] = v;
        T2[static_cast<std::size_t>(t)] = v * v;
    }

    // pattern: A = weight (1, or 0 at a wildcard), B = weight*value. both
    // reversed so the correlation over the sliding window is a plain
    // convolution. Csum = sum w*P^2 is the fixed, position-independent tail.
    std::vector<double> Arev(static_cast<std::size_t>(m), 0.0);
    std::vector<double> Brev(static_cast<std::size_t>(m), 0.0);
    double Csum = 0.0;
    for (int j = 0; j < m; ++j) {
        const char c = pat[static_cast<std::size_t>(j)];
        if (c == '?') continue;  // wildcard: weight 0, contributes nothing.
        const double v = static_cast<double>(c - 'a' + 1);
        const std::size_t rj = static_cast<std::size_t>(m - 1 - j);
        Arev[rj] = 1.0;
        Brev[rj] = v;
        Csum += v * v;
    }

    const std::vector<double> c1 = convolve(Arev, T2);  // sum w*T^2 over window
    const std::vector<double> c2 = convolve(Brev, T);   // sum w*P*T over window

    // M(i) = C1(i) - 2*C2(i) + Csum, read at index i+m-1. it's a nonnegative
    // integer bounded by 676*m; the float error is ~1e-6, so llround nails it.
    for (int i = 0; i + m <= n; ++i) {
        const std::size_t idx = static_cast<std::size_t>(i + m - 1);
        const double value = c1[idx] - 2.0 * c2[idx] + Csum;
        if (std::llround(value) == 0) res.push_back(i);
    }
    return res;
}

}  // namespace p0031
