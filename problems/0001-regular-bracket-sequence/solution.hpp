#pragma once

#include <cstdint>

namespace p0001 {

// the four length-2 blocks, by count:
//   c1 = "((",  c2 = "()",  c3 = ")(",  c4 = "))"
// glue all of them, in some order, into one string -- can it be a regular
// bracket sequence? 1 if yes, 0 if no.
//
// two facts settle it:
//   balance has to end at zero. every "((" is +2 and every "))" is -2, and the
//   neutral blocks cancel themselves -- so c1 == c4 or there's no hope.
//   a ")(" opens with a ')'. that ')' drives the running balance negative unless
//   it's already positive, and the only block that lifts balance ahead of
//   everything is "((". so any ")(" needs a "((" in front of it -- c1 > 0.
//   "()" is self-balanced. it never constrains a thing.
//
// order them "((" … ")(" … "()" … "))" and those two conditions are also enough:
// with c1 == c4 >= 1 the balance sits at 2*c1 through the middle, so every ")("
// dips one and recovers, and the trailing "))" walk it down to zero.
// c1 == c4 == 0 leaves only "()" and ")(" -- fine exactly when there's no ")(".
inline int regular(std::int64_t c1, std::int64_t c2, std::int64_t c3,
                   std::int64_t c4) {
    (void)c2;  // "()" is neutral -- any count of it is always placeable.
    return (c1 == c4 && (c3 == 0 || c1 > 0)) ? 1 : 0;
}

}  // namespace p0001
