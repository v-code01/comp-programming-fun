"""correctness gate for the digit DP -- differential vs the brute oracle, plus the
enumerated edges. runs under pytest (the test_* functions) and standalone (the
__main__ block prints "N/N checks passed" and exits 0 on a clean sweep, 1 on any
diff).

the brute is obviously right and dead slow. the DP is fast and has to prove it
agrees on thousands of small ranges -- if they ever disagree, the DP is wrong,
full stop. no mismatch is acceptable.
"""

import random
import sys

from solution import count, _count_upto
from reference import brute_count


def test_statement_examples():
    # 1..100, at most one distinct digit: 1..9 (nine singles) plus 11,22,…,99
    # (nine repdigits) plus 100? no -- 100 is {1,0}, two distinct. so 18.
    assert count(1, 100, 1) == 18
    # k=10 accepts everything -- only ten digit-values exist, nothing can exceed.
    assert count(1, 100, 10) == 100
    # single point.
    assert count(7, 7, 1) == 1
    # 100 alone has two distinct digits -- k=1 rejects it, k=2 keeps it.
    assert count(100, 100, 1) == 0
    assert count(100, 100, 2) == 1


def test_edges():
    # L == R across the whole space.
    for x in (1, 9, 10, 11, 99, 100, 101, 1000, 999999999999999999, 10 ** 18):
        for k in range(1, 11):
            assert count(x, x, k) == (1 if len(set(str(x))) <= k else 0)

    # K = 10 -> the entire range, always.
    assert count(1, 10 ** 6, 10) == 10 ** 6
    # K = 1 over 1..9 -> all nine single-digit numbers.
    assert count(1, 9, 1) == 9
    # ranges straddling powers of ten -- where digit-length changes, the classic
    # off-by-one trap for prefix counts.
    for lo, hi in [(1, 10), (9, 11), (99, 101), (999, 1001), (9999, 10001)]:
        for k in range(1, 11):
            assert count(lo, hi, k) == brute_count(lo, hi, k)

    # L = 1 lower boundary -- f(L-1) = f(0), and 0 must count as zero-distinct.
    assert _count_upto(0, 1) == 1  # just {0}, zero distinct <= 1.
    assert _count_upto(-1, 5) == 0  # nothing at or below -1.
    for hi in (5, 100, 12345):
        for k in range(1, 11):
            assert count(1, hi, k) == brute_count(1, hi, k)


def test_prefix_decomposition():
    # f(R) - f(L-1) must equal a direct brute over [L, R], wide but shallow so the
    # brute stays cheap. cross-checks the decomposition itself, not just count().
    rng = random.Random(20240702)
    for _ in range(300):
        a = rng.randint(1, 20000)
        b = rng.randint(1, 20000)
        L, R = min(a, b), max(a, b)
        k = rng.randint(1, 10)
        assert _count_upto(R, k) - _count_upto(L - 1, k) == brute_count(L, R, k)


def test_big_r_structural():
    # 1e18 sanity without a brute: k=10 is the full range width, and the count
    # can only rise as k rises (a looser bound never drops a number). monotone.
    R = 10 ** 18
    assert count(1, R, 10) == R  # every integer 1..1e18 has <= 10 distinct digits.
    prev = -1
    for k in range(1, 11):
        cur = count(1, R, k)
        assert cur >= prev  # monotone non-decreasing in k.
        assert 0 <= cur <= R
        prev = cur
    # repdigits 1..1e18 with k=1: 9 per length (d…d, d in 1..9) for lengths 1..18,
    # plus 1e18 itself which is {1,0}, not a repdigit. 9 * 18 = 162.
    assert count(1, R, 1) == 162


def _differential(rng, trials):
    """thousands of small random ranges -- DP vs brute. returns (checks, mismatches),
    printing every disagreement with its exact case and the two answers."""
    checks = 0
    mismatches = 0
    for _ in range(trials):
        L = rng.randint(1, 5000)
        span = rng.randint(0, 300)  # keep R - L small so the brute is cheap.
        R = L + span
        k = rng.randint(1, 10)
        got = count(L, R, k)
        want = brute_count(L, R, k)
        checks += 1
        if got != want:
            mismatches += 1
            print(f"MISMATCH  L={L} R={R} K={k}  dp={got}  brute={want}")
    return checks, mismatches


def test_differential():
    rng = random.Random(0xC0FFEE)
    checks, mismatches = _differential(rng, 4000)
    assert mismatches == 0, f"{mismatches} mismatches over {checks} checks"


if __name__ == "__main__":
    total = 0
    passed = 0
    fixed = [
        ("statement examples", test_statement_examples),
        ("edges", test_edges),
        ("prefix decomposition", test_prefix_decomposition),
        ("big R structural", test_big_r_structural),
    ]
    for name, fn in fixed:
        total += 1
        try:
            fn()
            passed += 1
        except AssertionError as e:
            print(f"FAIL [{name}]: {e}")

    # the differential sweep counts as one big check block -- report its internals
    # too so a green run shows the volume it actually cleared.
    rng = random.Random(0xC0FFEE)
    checks, mismatches = _differential(rng, 4000)
    total += 1
    if mismatches == 0:
        passed += 1
        print(f"differential: {checks}/{checks} random ranges agree, 0 diffs")
    else:
        print(f"differential: {mismatches} mismatches over {checks} checks")

    print(f"{passed}/{total} checks passed")
    sys.exit(0 if passed == total else 1)
