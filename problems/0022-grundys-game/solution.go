package main

import (
	"bufio"
	"fmt"
	"os"
)

// grundy's game. n heaps. a move picks one heap and splits it into two
// non-empty heaps of DIFFERENT sizes -- size s becomes (i, s-i) with
// 1 <= i < s-i. only heaps of size >= 3 can move at all: size 1 has nothing to
// split, and size 2 only splits into (1,1), which is equal, so it's dead too.
// last player to move wins -- normal play.
//
// sprague-grundy settles the whole thing. each heap is its own independent
// game. the value of the full position is the xor of the per-heap grundy
// values, and the first player wins iff that xor is nonzero.
//
// G(1)=G(2)=0 -- neither can move. for s >= 3:
//
//	G(s) = mex over 1<=i<s-i of { G(i) xor G(s-i) }
//
// here's the part people skip: grundy's-game values have no known closed form.
// they're only *conjectured* eventually periodic -- never proven, not in 100+
// years of people looking. so you can't table them or shortcut them. you
// compute every one, honestly, from the mex.
//
// the cost, stated plainly:
//   - reading input is Omega(n). you can't decide without seeing the heaps.
//   - building G(1..maxA) by the definition touches ~s/2 splits at each s, so
//     it's O(maxA^2) -- about 4.5e6 xor/mex ops at maxA=3000. that's the price
//     of no formula. no shortcut exists, so this is the floor, not sloppiness.
//   - after that the verdict is one pass: O(n) xors.
//
// total: O(maxA^2 + n). the square term dominates.

// grundy computes G(0..maxA). g[0] is unused; g[1]=g[2]=0 fall out of the
// zero value. everything past that is the mex of the split values.
func grundy(maxA int) []int {
	g := make([]int, maxA+1)

	// one seen-buffer, reused across every s -- no per-s allocation churn.
	// mex(s) can't exceed the number of split options, which is (s-1)/2 < s,
	// so a buffer this wide indexes safely for every s.
	seen := make([]bool, maxA+2)

	for s := 3; s <= maxA; s++ {
		half := (s - 1) / 2 // i runs 1..half -- that many split options.

		for i := 1; i <= half; i++ {
			v := g[i] ^ g[s-i]
			// mex tops out at half (the set can be at most {0..half-1..half}).
			// any split value bigger than half can't be the mex, so it can't
			// block it either -- drop it and keep the wipe cheap.
			if v <= half {
				seen[v] = true
			}
		}

		mex := 0
		for seen[mex] {
			mex++
		}
		g[s] = mex

		// wipe only what we touched -- 0..half. keeps each s at O(s), not
		// O(maxA), so the whole build stays O(maxA^2).
		for i := 0; i <= half; i++ {
			seen[i] = false
		}
	}
	return g
}

// firstWins is the solver. xor the grundy value of every heap; the first
// player wins iff the result is nonzero. this is the sprague-grundy verdict,
// nothing more.
func firstWins(heaps []int) bool {
	maxA := 1
	for _, h := range heaps {
		if h > maxA {
			maxA = h
		}
	}

	g := grundy(maxA)

	x := 0
	for _, h := range heaps {
		x ^= g[h]
	}
	return x != 0
}

func main() {
	reader := bufio.NewReader(os.Stdin)

	var n int
	if _, err := fmt.Fscan(reader, &n); err != nil {
		return // empty input -- nothing to decide.
	}

	heaps := make([]int, n)
	for i := 0; i < n; i++ {
		fmt.Fscan(reader, &heaps[i])
	}

	if firstWins(heaps) {
		fmt.Println("First")
	} else {
		fmt.Println("Second")
	}
}
