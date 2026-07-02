#include <cstdint>
#include <cstdio>

#include "common/testing.hpp"
#include "reference.hpp"
#include "solution.hpp"

int main() {
    using namespace p0001;

    // the three examples from the statement.
    CPF_EQ(regular(3, 1, 4, 3), 1);
    CPF_EQ(regular(0, 0, 0, 0), 1);
    CPF_EQ(regular(1, 2, 3, 4), 0);

    // edges the formula has to get right:
    CPF_EQ(regular(0, 0, 1, 0), 0);   // a lone ")(" starts negative.
    CPF_EQ(regular(0, 5, 1, 0), 0);   // "()" can't rescue a leading ")(".
    CPF_EQ(regular(1, 0, 1, 1), 1);   // one "((" covers the ")(".
    CPF_EQ(regular(2, 0, 0, 1), 0);   // unbalanced -- c1 != c4.
    CPF_EQ(regular(5, 0, 0, 5), 1);   // balanced, no ")(" to worry about.

    // full 1e9 counts -- the check is O(1) and never adds, so no overflow.
    CPF_EQ(regular(1000000000LL, 0, 1000000000LL, 1000000000LL), 1);
    CPF_EQ(regular(1000000000LL, 0, 0, 999999999LL), 0);

    // the whole branch space lives in a handful of small values, so don't
    // sample it -- exhaust it. every combination of counts in 0..6 has to match
    // the search that actually tries to build the string.
    ReferenceSolver ref;
    int cases = 0, diffs = 0, first_a = -1, first_b = -1, first_c = -1,
        first_d = -1;
    for (int a = 0; a <= 6; ++a)
        for (int b = 0; b <= 6; ++b)
            for (int c = 0; c <= 6; ++c)
                for (int d = 0; d <= 6; ++d) {
                    ++cases;
                    if (regular(a, b, c, d) != ref.solve(a, b, c, d)) {
                        if (diffs == 0) {
                            first_a = a;
                            first_b = b;
                            first_c = c;
                            first_d = d;
                        }
                        ++diffs;
                    }
                }
    CPF_CHECK(diffs == 0);
    if (diffs) {
        std::printf("  first diff at (%d,%d,%d,%d)\n", first_a, first_b, first_c,
                    first_d);
    }
    std::printf("exhaustive grid 0..6^4: %d cases, %d diffs\n", cases, diffs);

    return cpf::report();
}
