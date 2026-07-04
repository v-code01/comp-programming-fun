package main

import (
	"sort"
	"strconv"
	"strings"
)

// bruteFirstWins -- the independent ground truth. no grundy, no xor, no theory
// borrowed from the solver. just the raw question: can the player to move force
// a win?
//
// a position is the whole multiset of heaps. a move splits one heap of size h
// into two unequal non-empty parts (i, h-i), 1 <= i < h-i. no move available
// means you just lost -- normal play. a position is WINNING iff some move hands
// the opponent a LOSING position. that's minimax, straight.
//
// memoized on the sorted-multiset key so identical positions aren't re-solved,
// but the state space still explodes -- this is exponential. feed it tiny
// positions only. it exists to check the solver, not to replace it.
func bruteFirstWins(heaps []int) bool {
	h := append([]int(nil), heaps...)
	sort.Ints(h)
	return winning(h, map[string]bool{})
}

// key canonicalizes a position: sort, then join. two positions with the same
// heap multiset share a key -- that's what makes the memo hit.
func key(heaps []int) string {
	parts := make([]string, len(heaps))
	for i, v := range heaps {
		parts[i] = strconv.Itoa(v)
	}
	return strings.Join(parts, ",")
}

// winning returns true iff the player to move on this sorted position wins.
func winning(heaps []int, memo map[string]bool) bool {
	k := key(heaps)
	if v, ok := memo[k]; ok {
		return v
	}

	res := false
	for idx, h := range heaps {
		if h < 3 {
			continue // nothing splits into two unequal non-empty parts.
		}
		for i := 1; i < h-i; i++ { // i < h-i forces the two parts unequal.
			next := make([]int, 0, len(heaps)+1)
			for j, v := range heaps {
				if j != idx {
					next = append(next, v)
				}
			}
			next = append(next, i, h-i)
			sort.Ints(next)

			// hand the opponent this position. if it's a loss for them, we win.
			if !winning(next, memo) {
				res = true
				break
			}
		}
		if res {
			break
		}
	}

	memo[k] = res
	return res
}
