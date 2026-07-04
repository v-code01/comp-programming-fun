package main

// diameter of a point set -- the max squared distance between any two points.
//
// the naive move is to look at every pair. that's O(n^2) -- 1e10 compares at
// n=1e5, dead on arrival. the fact that saves you: the farthest pair both sit
// on the convex hull. nothing interior can be an endpoint of the diameter --
// push any interior point outward along the segment to the far point and the
// distance only grows, so the extreme already lives on the boundary.
//
// so the plan is two moves:
//   build the hull -- Andrew's monotone chain, O(n log n). sort by (x, y),
//   sweep the lower chain then the upper, popping every right turn and every
//   collinear point. what's left is the strictly-convex boundary, CCW.
//   walk the hull with rotating calipers -- O(h). for each edge, slide the
//   antipodal vertex forward to the one farthest from that edge, and the
//   diameter shows up as a distance to one of those antipodal vertices.
//
// LOWER BOUND -- what's forced, honestly:
//   reading n points is Omega(n). you can't answer without looking at them.
//   the hull needs a sort, and sorting is Omega(n log n) in the comparison
//   model. the calipers walk is O(h) <= O(n) on top -- free next to the sort.
//   so O(n log n) total, and that's the floor here, not a loose bound. the
//   farthest-pair-on-hull fact is the whole reason we dodge the O(n^2) scan --
//   without it you're back to all pairs.
//
// everything is integer. coords fit in 1e6, so a coordinate difference fits in
// 2e6, a cross product in 2*(2e6)^2 = 8e12, a squared distance in 8e12 -- all
// well inside int64. no floats touch the answer. no floats touch the hull.

import (
	"bufio"
	"os"
	"sort"
	"strconv"
)

// Point is one integer lattice point. it's comparable, so == dedups directly.
type Point struct {
	x, y int64
}

// area2 is twice the signed area of triangle (a, b, c) -- the cross product of
// (b-a) and (c-a). positive means c is left of the ray a->b (a CCW turn).
// this is the only orientation primitive the hull and calipers need.
func area2(a, b, c Point) int64 {
	return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x)
}

// distSq is the squared Euclidean distance -- exact, no sqrt, no float. this is
// the answer's currency.
func distSq(a, b Point) int64 {
	dx := a.x - b.x
	dy := a.y - b.y
	return dx*dx + dy*dy
}

func maxi64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

// convexHull returns the strictly-convex boundary in CCW order. collinear
// points on an edge are dropped -- a diameter endpoint is always a corner, so
// the interior of an edge never matters.
//
// it copies first -- the caller's slice keeps its order, which the brute oracle
// leans on during differential tests.
//
// degeneracies fall out for free. duplicates die in the dedup pass. one unique
// point returns one point. two unique points return the pair. all-collinear
// collapses to the two extremes -- the sweep pops every point between them.
func convexHull(in []Point) []Point {
	pts := make([]Point, len(in))
	copy(pts, in)

	sort.Slice(pts, func(i, j int) bool {
		if pts[i].x != pts[j].x {
			return pts[i].x < pts[j].x
		}
		return pts[i].y < pts[j].y
	})

	// dedup in place -- duplicate points would wreck the collinearity test and
	// they can never be a distinct diameter endpoint anyway.
	m := 0
	for i := range pts {
		if m == 0 || pts[m-1] != pts[i] {
			pts[m] = pts[i]
			m++
		}
	}
	pts = pts[:m]

	n := len(pts)
	if n <= 2 {
		// 0, 1, or 2 unique points -- the hull is just the points. calipers
		// reads it straight.
		return pts
	}

	hull := make([]Point, 0, 2*n)

	// lower chain, left to right.
	for i := 0; i < n; i++ {
		for len(hull) >= 2 && area2(hull[len(hull)-2], hull[len(hull)-1], pts[i]) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, pts[i])
	}

	// upper chain, right to left. the +1 floor keeps the lower chain's last
	// point pinned while the upper builds on top of it.
	lower := len(hull) + 1
	for i := n - 2; i >= 0; i-- {
		for len(hull) >= lower && area2(hull[len(hull)-2], hull[len(hull)-1], pts[i]) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, pts[i])
	}

	// the last point is the start point repeated -- drop it.
	return hull[:len(hull)-1]
}

// diameterSq walks the hull with rotating calipers and returns the max squared
// distance between any two vertices -- which, since the hull holds the extremes,
// is the diameter of the whole set.
//
// the loop is bug-prone exactly where the hull is small. a 1-vertex hull (all
// points identical) has diameter 0. a 2-vertex hull (n=2, or everything
// collinear) is one segment -- answer is that segment's length. both are
// handled before the calipers ever spins, so it never indexes a phantom edge
// and never asks for the area of a degenerate triangle.
//
// for h >= 3: for each edge (i, i+1) advance the antipodal vertex j while the
// next vertex sits farther from the edge's line -- area2 is that distance times
// a constant, so the compare is exact integer, no division, no zero divide. j
// only moves forward and laps the hull a bounded number of times, so the whole
// walk is O(h). checking both edge endpoints against j at every step is what
// keeps parallel supporting edges from hiding the true pair.
func diameterSq(h []Point) int64 {
	n := len(h)
	if n <= 1 {
		return 0
	}
	if n == 2 {
		return distSq(h[0], h[1])
	}

	best := int64(0)
	j := 1
	for i := 0; i < n; i++ {
		inext := (i + 1) % n
		// slide j to the vertex farthest from edge (i, inext).
		for area2(h[i], h[inext], h[(j+1)%n]) > area2(h[i], h[inext], h[j]) {
			j = (j + 1) % n
		}
		best = maxi64(best, distSq(h[i], h[j]))
		best = maxi64(best, distSq(h[inext], h[j]))
	}
	return best
}

// solveDiameterSq is the whole pipeline -- hull, then calipers.
func solveDiameterSq(pts []Point) int64 {
	return diameterSq(convexHull(pts))
}

// readInt pulls one signed integer off the reader, skipping leading whitespace.
// hand-rolled -- fmt.Fscan at n=1e5 is the slow path we don't need.
func readInt(r *bufio.Reader) int64 {
	c, err := r.ReadByte()
	for err == nil && (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
		c, err = r.ReadByte()
	}
	neg := false
	if c == '-' {
		neg = true
		c, err = r.ReadByte()
	}
	var v int64
	for err == nil && c >= '0' && c <= '9' {
		v = v*10 + int64(c-'0')
		c, err = r.ReadByte()
	}
	if neg {
		return -v
	}
	return v
}

func main() {
	r := bufio.NewReaderSize(os.Stdin, 1<<20)
	n := int(readInt(r))
	pts := make([]Point, n)
	for i := 0; i < n; i++ {
		pts[i].x = readInt(r)
		pts[i].y = readInt(r)
	}

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	w.WriteString(strconv.FormatInt(solveDiameterSq(pts), 10))
	w.WriteByte('\n')
}
