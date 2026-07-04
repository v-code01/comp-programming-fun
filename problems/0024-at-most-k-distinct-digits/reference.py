"""brute oracle -- walk every x in [L, R] and count its distinct digits directly.

no cleverness on purpose. set(str(x)) is the distinct-digit set of x written with
no leading zeros -- str never pads, so 100 is "100" -> {'1','0','0'} -> {'1','0'},
two distinct. this is the ground truth the digit DP has to match. O(R - L + 1),
only ever called on tiny ranges in the tests.
"""


def brute_count(L: int, R: int, K: int) -> int:
    """count x in [L, R] with at most K distinct decimal digits, the slow honest way."""
    total = 0
    for x in range(L, R + 1):
        if len(set(str(x))) <= K:
            total += 1
    return total


if __name__ == "__main__":
    import sys

    L, R, K = (int(t) for t in sys.stdin.read().split())
    print(brute_count(L, R, K))
