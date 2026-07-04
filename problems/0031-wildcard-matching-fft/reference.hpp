#pragma once

#include <string>
#include <vector>

namespace p0031 {

// the oracle. independent of the FFT, and obviously right: for every start i,
// walk all m pattern chars and demand each is a '?' or equals the text char
// under it. O(n*m). small inputs only -- this is the ground truth the FFT has
// to reproduce on thousands of random cases, not a solver you'd ship.
//
// returns the matching starts ascending, same contract as matches().
inline std::vector<int> reference(const std::string& text,
                                  const std::string& pat) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pat.size());
    std::vector<int> res;
    if (m == 0 || m > n) return res;

    for (int i = 0; i + m <= n; ++i) {
        bool ok = true;
        for (int j = 0; j < m; ++j) {
            const char pc = pat[static_cast<std::size_t>(j)];
            if (pc == '?') continue;  // wildcard matches anything.
            if (pc != text[static_cast<std::size_t>(i + j)]) {
                ok = false;
                break;
            }
        }
        if (ok) res.push_back(i);
    }
    return res;
}

}  // namespace p0031
