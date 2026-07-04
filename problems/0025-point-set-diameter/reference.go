package main

// the oracle. no hull, no calipers, no cleverness -- just look at every pair and
// keep the biggest squared distance. O(n^2). it's obviously right because it's
// the literal definition of the diameter, and being obviously right is its whole
// job. the differential test trusts this and doubts the fast path.
//
// only ever run on small sets -- at n=1e5 this is 1e10 compares and you'd wait
// all day. that's fine. it exists to catch calipers bugs on n<=40, where the
// degenerate hulls that break rotating calipers actually live.
func bruteDiameterSq(pts []Point) int64 {
	best := int64(0)
	for i := 0; i < len(pts); i++ {
		for j := i + 1; j < len(pts); j++ {
			if d := distSq(pts[i], pts[j]); d > best {
				best = d
			}
		}
	}
	return best
}
