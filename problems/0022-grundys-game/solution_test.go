package main

import (
	"math/rand/v2"
	"testing"
)

// TestExamples nails the statement cases and the hand-checked edges.
func TestExamples(t *testing.T) {
	cases := []struct {
		name  string
		heaps []int
		want  bool // true = "First" wins.
	}{
		// single unsplittable heaps -- no move exists, so the first player is
		// stuck. second wins both.
		{"single-1", []int{1}, false},
		{"single-2", []int{2}, false},

		// size 3 splits into (1,2), leaving the opponent stuck. first wins.
		{"single-3", []int{3}, true},

		// {1,3}: split the 3 into (1,2) -> {1,1,2}, all dead. opponent can't
		// move. first wins. G(1)^G(3) = 0^1 = 1 != 0 agrees.
		{"one-and-three", []int{1, 3}, true},

		// every heap dead -- first player can't move at all. second wins.
		{"all-twos", []int{2, 2, 2, 2}, false},
		{"ones-and-twos", []int{1, 1, 2, 2, 1}, false},

		// two size-3 heaps. G(3)^G(3) = 1^1 = 0 -> second. mirror strategy.
		{"two-threes", []int{3, 3}, false},

		// one big heap. G(5) = 2 != 0 -> first.
		{"single-5", []int{5}, true},

		// mixed multi-heap, cross-checked against the oracle below anyway.
		{"mixed", []int{4, 5, 6}, bruteFirstWins([]int{4, 5, 6})},
	}

	for _, c := range cases {
		if got := firstWins(c.heaps); got != c.want {
			t.Fatalf("%s: firstWins(%v) = %v, want %v", c.name, c.heaps, got, c.want)
		}
	}
}

// TestGrundyValues pins the first few grundy numbers computed by hand, so a
// regression in the mex loop can't hide behind the xor.
func TestGrundyValues(t *testing.T) {
	// G: index 0..12. hand-derived from the mex definition.
	//   G1=0 G2=0 G3=1 G4=0 G5=2 G6=1 G7=0 G8=2 G9=1 G10=0 G11=2 G12=1
	want := []int{0, 0, 0, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1}
	g := grundy(12)
	for s := 1; s <= 12; s++ {
		if g[s] != want[s] {
			t.Fatalf("G(%d) = %d, want %d", s, g[s], want[s])
		}
	}
}

// TestDifferential is the correctness gate. thousands of small random
// positions, solver vs the independent minimax oracle. any single mismatch
// prints the exact position and stops the world.
func TestDifferential(t *testing.T) {
	rng := rand.New(rand.NewPCG(0x9E3779B97F4A7C15, 0xD1B54A32D192ED03))

	const trials = 5000
	for iter := 0; iter < trials; iter++ {
		// n <= 4, a_i in [1,14] keeps the exponential oracle cheap while still
		// exercising real splits (heaps up to 14 branch plenty).
		n := 1 + rng.IntN(4)
		heaps := make([]int, n)
		for i := range heaps {
			heaps[i] = 1 + rng.IntN(14)
		}

		got := firstWins(heaps)
		want := bruteFirstWins(heaps)
		if got != want {
			t.Fatalf("mismatch on %v: solver=%v oracle=%v", heaps, got, want)
		}
	}
}

// BenchmarkSolve measures the real envelope: maxA=3000, n=1000. the O(maxA^2)
// grundy build dominates -- that's the number to watch.
func BenchmarkSolve(b *testing.B) {
	heaps := make([]int, 1000)
	for i := range heaps {
		heaps[i] = 1 + (i*7)%3000 // spread across [1,3000]; index 0 forces 3000.
	}
	heaps[0] = 3000 // pin maxA at the constraint ceiling.

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if firstWins(heaps) {
			b.SetBytes(1) // touch the result so the call can't be optimized away.
		}
	}
}
