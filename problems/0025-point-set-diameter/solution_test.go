package main

import (
	"math"
	"math/rand/v2"
	"testing"
)

// fmtPts renders a point set so a failing case is copy-pasteable back into a
// test. a diff you can't reproduce is a diff you can't kill.
func fmtPts(pts []Point) string {
	s := "["
	for i, p := range pts {
		if i > 0 {
			s += " "
		}
		s += "(" + itoa(p.x) + "," + itoa(p.y) + ")"
	}
	return s + "]"
}

func itoa(v int64) string {
	if v == 0 {
		return "0"
	}
	neg := v < 0
	if neg {
		v = -v
	}
	var b [24]byte
	i := len(b)
	for v > 0 {
		i--
		b[i] = byte('0' + v%10)
		v /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}

// the statement's worked examples, plus the degenerate corners spelled out.
func TestExamples(t *testing.T) {
	cases := []struct {
		name string
		pts  []Point
		want int64
	}{
		{"right-triangle-3-4-5", []Point{{0, 0}, {3, 0}, {0, 4}}, 25},
		{"unit-square", []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 2},
		{"collinear-span", []Point{{0, 0}, {2, 0}, {5, 0}, {1, 0}}, 25},
		{"n-equals-2", []Point{{1, 1}, {4, 5}}, 25},
		{"all-identical", []Point{{7, 7}, {7, 7}, {7, 7}}, 0},
		{"two-identical", []Point{{3, 3}, {3, 3}}, 0},
		{"vertical-line", []Point{{2, -3}, {2, 10}, {2, 4}}, 169},
		{"duplicates-plus-far", []Point{{0, 0}, {0, 0}, {0, 0}, {6, 8}}, 100},
		{"negatives", []Point{{-1000000, -1000000}, {1000000, 1000000}}, 8000000000000},
	}
	for _, c := range cases {
		if got := solveDiameterSq(c.pts); got != c.want {
			t.Fatalf("%s: got %d, want %d for %s", c.name, got, c.want, fmtPts(c.pts))
		}
	}
}

// the real gate. thousands of random sets over tight coordinate ranges -- small
// ranges are on purpose, they manufacture collinear triples and duplicates by
// the fistful, which is exactly where a rotating-calipers bug hides. every set
// gets checked against the brute oracle.
func TestDifferential(t *testing.T) {
	rng := rand.New(rand.NewPCG(0xC0FFEE, 0x1234567))

	check := func(tag string, pts []Point) {
		want := bruteDiameterSq(pts)
		got := solveDiameterSq(pts)
		if got != want {
			t.Fatalf("%s: got %d, want %d for %s", tag, got, want, fmtPts(pts))
		}
	}

	// tight random cloud -- ranges from 1 to 6 so points pile up and line up.
	for iter := 0; iter < 6000; iter++ {
		n := 2 + rng.IntN(39) // 2..40
		span := int64(1 + rng.IntN(6))
		pts := make([]Point, n)
		for i := range pts {
			pts[i] = Point{rng.Int64N(2*span+1) - span, rng.Int64N(2*span+1) - span}
		}
		check("random-tight", pts)
	}

	// wider random cloud -- bigger coords, fewer collisions, real hulls.
	for iter := 0; iter < 3000; iter++ {
		n := 2 + rng.IntN(39)
		pts := make([]Point, n)
		for i := range pts {
			pts[i] = Point{rng.Int64N(200001) - 100000, rng.Int64N(200001) - 100000}
		}
		check("random-wide", pts)
	}

	// all-collinear -- points scattered on one random line. the hull collapses
	// to a segment and calipers must not spin on it.
	for iter := 0; iter < 2000; iter++ {
		n := 2 + rng.IntN(39)
		ax := rng.Int64N(21) - 10
		ay := rng.Int64N(21) - 10
		bx := rng.Int64N(2001) - 1000
		by := rng.Int64N(2001) - 1000
		pts := make([]Point, n)
		for i := range pts {
			t := rng.Int64N(41) - 20
			pts[i] = Point{bx + t*ax, by + t*ay}
		}
		check("collinear", pts)
	}

	// on a circle -- integer-rounded angles give a big, near-regular hull with
	// nearly every point a vertex. this is where the antipodal walk earns its
	// keep, and where parallel-ish supporting lines show up.
	for iter := 0; iter < 2000; iter++ {
		n := 3 + rng.IntN(30)
		r := float64(50 + rng.IntN(950))
		cx := int64(rng.IntN(201) - 100)
		cy := int64(rng.IntN(201) - 100)
		pts := make([]Point, n)
		for i := range pts {
			a := 2 * math.Pi * float64(i) / float64(n)
			pts[i] = Point{cx + int64(math.Round(r*math.Cos(a))), cy + int64(math.Round(r*math.Sin(a)))}
		}
		check("on-circle", pts)
	}

	// duplicate-heavy -- draw from a tiny alphabet of points so most are repeats.
	for iter := 0; iter < 2000; iter++ {
		n := 2 + rng.IntN(39)
		alpha := make([]Point, 1+rng.IntN(4)) // 1..4 distinct points
		for i := range alpha {
			alpha[i] = Point{rng.Int64N(2001) - 1000, rng.Int64N(2001) - 1000}
		}
		pts := make([]Point, n)
		for i := range pts {
			pts[i] = alpha[rng.IntN(len(alpha))]
		}
		check("duplicate-heavy", pts)
	}

	// regular polygons -- even side counts give exactly parallel opposite edges,
	// the classic trap for the "advance while strictly greater" caliper step.
	for _, sides := range []int{3, 4, 5, 6, 8, 10, 12, 20, 50, 100} {
		r := 1000.0
		pts := make([]Point, sides)
		for i := range pts {
			a := 2 * math.Pi * float64(i) / float64(sides)
			pts[i] = Point{int64(math.Round(r * math.Cos(a))), int64(math.Round(r * math.Sin(a)))}
		}
		check("regular-polygon", pts)
	}
}

// n=1e5 on a large circle -- the worst case. the hull is enormous (thousands of
// vertices), so this measures the calipers walk, not a hull that shrank to a
// handful of corners.
func BenchmarkSolve(b *testing.B) {
	const n = 100000
	r := 999999.0
	pts := make([]Point, n)
	for i := range pts {
		a := 2 * math.Pi * float64(i) / float64(n)
		pts[i] = Point{int64(math.Round(r * math.Cos(a))), int64(math.Round(r * math.Sin(a)))}
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = solveDiameterSq(pts)
	}
}
