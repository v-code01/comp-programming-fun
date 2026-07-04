package main

import "sort"

// the oracle. dumb on purpose -- it's the thing the fast solution has to agree
// with. it enumerates all 2^n subsets, XORs each, dedups, sorts. exponential, so
// n stays tiny (<= 16). if the basis ever disagrees with this, the basis is wrong.

// bruteDistinct returns every distinct achievable XOR of a, sorted ascending.
// O(2^n * n). empty subset contributes 0 at mask == 0.
func bruteDistinct(a []uint64) []uint64 {
	n := len(a)
	seen := make(map[uint64]struct{})
	total := 1 << uint(n)
	for mask := 0; mask < total; mask++ {
		var x uint64
		for i := 0; i < n; i++ {
			if mask&(1<<uint(i)) != 0 {
				x ^= a[i]
			}
		}
		seen[x] = struct{}{}
	}
	out := make([]uint64, 0, len(seen))
	for v := range seen {
		out = append(out, v)
	}
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}

// bruteKth returns the k-th smallest distinct achievable XOR (k is 1-indexed), or
// -1 when k runs past the count. the reference answer the fast path is graded on.
func bruteKth(a []uint64, k uint64) int64 {
	d := bruteDistinct(a)
	if k == 0 || k > uint64(len(d)) {
		return -1
	}
	return int64(d[k-1])
}
