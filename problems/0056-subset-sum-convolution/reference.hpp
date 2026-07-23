#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace p0056 {

// dumb and obviously correct, small k only. this is the definition, typed out:
// c[S] = sum over T subset of S of a[T]*b[S\T], mod p. no zeta, no ranks, no
// covering product -- just walk the submasks of S with the classic
//   for (T = S; ; T = (T-1) & S)
// enumeration and accumulate a[T]*b[S^T]. S\T is S with T's bits cleared, which
// is S ^ T exactly when T is a submask of S -- and this loop only ever visits
// submasks. it is O(3^k) over all S (each bit lands in T, in S\T, or outside S),
// so it stays cheap only while k is small; 3^12 ~ 5.3e5 is the practical ceiling.
//
// the whole point of this file is independence: it shares no line of machinery
// with solve(), so agreeing on thousands of random instances is real evidence,
// not two copies of the same bug nodding at each other.
inline std::vector<std::uint32_t> reference(const std::vector<std::uint32_t>& a,
                                            const std::vector<std::uint32_t>& b,
                                            int k) {
    const std::size_t size = static_cast<std::size_t>(1) << k;
    std::vector<std::uint32_t> c(size);
    for (std::size_t S = 0; S < size; ++S) {
        std::uint64_t acc = 0;
        // submask walk. terminates: T strictly decreases each step within the
        // bits of S, and the T==0 break fires after the empty submask is counted.
        for (std::size_t T = S;; T = (T - 1) & S) {
            acc += static_cast<std::uint64_t>(a[T]) * b[S ^ T] % 998244353ULL;
            if (T == 0) break;
        }
        c[S] = static_cast<std::uint32_t>(acc % 998244353ULL);
    }
    return c;
}

}  // namespace p0056
