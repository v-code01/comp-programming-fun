"""at most k distinct digits -- count x in [L, R] whose decimal spelling uses
at most K different digit-values.

the honest way is to walk every x from L to R and count its distinct digits.
that's R - L + 1 steps, and R goes to 1e18 -- dead. so we don't count numbers,
we count digit-strings. fix an upper bound N and ask: how many x in [0, N] spell
out at most K distinct digits? call that f(N). then

    answer = f(R) - f(L-1).

both f(R) and f(L-1) count from 0, and subtracting drops the [0, L-1] tail --
standard prefix-count trick, the whole range [L, R] minus the range below L.

the DP walks N one digit at a time, most-significant first, carrying four things:
  pos     -- which digit index we're placing.
  mask    -- 10-bit set of digit-values used so far (bit d means d appeared).
  tight   -- are we still hugging N's own prefix? if yes, this digit caps at
             N[pos]; once we place something strictly smaller, tight drops and
             every lower position is free (0..9).
  started -- have we placed a non-zero digit yet? this is the whole leading-zero
             story, see below.

why 'started' -- the leading-zero trap:
  100 spells {1, 0} -- two distinct digits. but 007 is just 7, one digit. the
  leading zeros in front of the first real digit are padding, not the digit '0'.
  if we blindly OR bit 0 for every zero we place, then 100 walks 1, then 0, then
  0 and reports {1, 0}… which is right, but 7 padded to width 3 as 007 would
  report {0, 7} -- wrong, that '0' never got written.
  the fix: while started is false we are still in the padding. a 0 placed there
  sets nothing -- mask stays empty, we stay not-started. the first non-zero digit
  flips started to true and from then on every digit, including any 0, is real
  and gets OR'd into mask. so 100 = [1 flips started, mask={1}][0 real, mask={1,0}]
  [0 real, mask={1,0}] -> {1,0}, two distinct. and 7 never pads because we only
  descend as many positions as N actually has -- pos runs over N's own length.
  the number 0 itself finishes with started still false and mask empty; we count
  it as valid (zero distinct <= K), which is what f(N) wants for the [0, N] range.

accept rule at the end: popcount(mask) <= K. an all-zero mask (the number 0, or
the not-yet-started path) has popcount 0 <= K -- always fine.

complexity. the state is (pos, mask, tight, started): len * 2^10 * 2 * 2. len is
the digit count of N, <= 19 for 1e18. transitions try <= 10 digits each. so each
f(N) is O(len * 2^10 * 10) memo fills -- a few hundred thousand, constant in the
*value* of N because it only depends on N's length and the fixed 10-digit alphabet.
two calls for the answer. wall time is microseconds.

lower bound, honestly. you have to read L, R, K -- that's the input, unavoidable.
you have to look at N's ~19 digits to know what you're bounding under; you can't
count a range without knowing its top. and the mask lives in a fixed 2^10 space
because there are ten decimal digits, not a variable you can shrink. so the work
is Theta(len + 2^10) per bound at the floor and our O(len * 2^10 * 10) is that
same fixed budget times a small constant -- there's no sublinear dodge, the digits
and the alphabet are both irreducible.
"""

from functools import lru_cache


def _count_upto(n: int, k: int) -> int:
    """f(n) -- how many x in [0, n] spell at most k distinct digits.

    n >= 0. returns the count. digit DP, memoized on (pos, mask, tight, started).
    """
    if n < 0:
        return 0  # empty prefix -- nothing at or below a negative bound.

    digits = [int(c) for c in str(n)]  # most-significant first.
    length = len(digits)

    @lru_cache(maxsize=None)
    def walk(pos: int, mask: int, tight: bool, started: bool) -> int:
        # placed every position -- this path spelled one whole number. accept it
        # iff its distinct-digit count fits under k. bin(mask).count('1') is the
        # popcount; mask==0 (the number 0 or a still-padding path) gives 0 <= k.
        if pos == length:
            return 1 if bin(mask).count("1") <= k else 0

        # tight caps this digit at n's own digit here -- go past it and we'd
        # exceed n. loose means the whole 0..9 alphabet is open.
        hi = digits[pos] if tight else 9
        total = 0
        for d in range(hi + 1):
            next_tight = tight and (d == hi)  # still hugging n only on the cap.
            if not started and d == 0:
                # padding zero -- not a real digit. mask and started stay put.
                total += walk(pos + 1, mask, next_tight, False)
            else:
                # first non-zero flips started; every real digit joins the mask.
                total += walk(pos + 1, mask | (1 << d), next_tight, True)
        return total

    result = walk(0, 0, True, False)
    walk.cache_clear()  # memo is keyed to *this* n -- don't leak into the next call.
    return result


def count(L: int, R: int, K: int) -> int:
    """count x in [L, R] with at most K distinct decimal digits.

    prefix subtraction: f(R) minus f(L-1) strips the [0, L-1] tail and leaves
    exactly [L, R]. K >= 10 accepts everything (only ten digits exist), K = 1
    keeps single-digit numbers and repdigits like 777.
    """
    return _count_upto(R, K) - _count_upto(L - 1, K)


if __name__ == "__main__":
    import sys

    L, R, K = (int(t) for t in sys.stdin.read().split())
    print(count(L, R, K))
