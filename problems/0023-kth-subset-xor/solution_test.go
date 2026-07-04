package main

import (
	"math/rand/v2"
	"testing"
)

// TestExamples pins the statement's own numbers and the sharp edges around them.
func TestExamples(t *testing.T) {
	// [1,2] spans {0,1,2,3}. four values, then -1.
	b := newBasis([]uint64{1, 2})
	want := []int64{0, 1, 2, 3, -1}
	for i, w := range want {
		k := uint64(i + 1)
		if got := b.kth(k); got != w {
			t.Fatalf("[1,2] k=%d: got %d want %d", k, got, w)
		}
	}

	// all zeros -- the only reachable value is 0. rank 0, count 1.
	b = newBasis([]uint64{0, 0, 0})
	if got := b.kth(1); got != 0 {
		t.Fatalf("zeros k=1: got %d want 0", got)
	}
	if got := b.kth(2); got != -1 {
		t.Fatalf("zeros k=2: got %d want -1", got)
	}

	// empty array -- still {0}. the empty subset is always there.
	b = newBasis(nil)
	if got := b.kth(1); got != 0 {
		t.Fatalf("empty k=1: got %d want 0", got)
	}
	if got := b.kth(2); got != -1 {
		t.Fatalf("empty k=2: got %d want -1", got)
	}

	// linearly dependent set: 3 = 1 XOR 2, so rank is 2, not 3. still {0,1,2,3}.
	b = newBasis([]uint64{1, 2, 3})
	want = []int64{0, 1, 2, 3, -1}
	for i, w := range want {
		k := uint64(i + 1)
		if got := b.kth(k); got != w {
			t.Fatalf("[1,2,3] k=%d: got %d want %d", k, got, w)
		}
	}

	// duplicates collapse -- [5,5,5] has rank 1, values {0,5}.
	b = newBasis([]uint64{5, 5, 5})
	if got := b.kth(1); got != 0 {
		t.Fatalf("[5,5,5] k=1: got %d want 0", got)
	}
	if got := b.kth(2); got != 5 {
		t.Fatalf("[5,5,5] k=2: got %d want 5", got)
	}
	if got := b.kth(3); got != -1 {
		t.Fatalf("[5,5,5] k=3: got %d want -1", got)
	}

	// full rank on a wide bit -- a value near 2^60 still sorts and indexes right.
	hi := uint64(1) << 59
	b = newBasis([]uint64{hi, 1})
	// {0, 1, hi, hi+1} ascending.
	wantWide := []int64{0, 1, int64(hi), int64(hi + 1), -1}
	for i, w := range wantWide {
		k := uint64(i + 1)
		if got := b.kth(k); got != w {
			t.Fatalf("wide k=%d: got %d want %d", k, got, w)
		}
	}
}

// TestDifferential is the correctness gate. thousands of random small arrays,
// every valid k plus a few past the end -- basis must match the 2^n oracle exactly.
func TestDifferential(t *testing.T) {
	rng := rand.New(rand.NewPCG(0xC0FFEE, 0x1234567))
	const trials = 4000
	for tr := 0; tr < trials; tr++ {
		n := rng.IntN(13)           // 0..12 elements
		valBits := 1 + rng.IntN(12) // values < 2^valBits, up to 2^12
		mask := (uint64(1) << uint(valBits)) - 1
		a := make([]uint64, n)
		for i := range a {
			a[i] = rng.Uint64() & mask
		}

		d := bruteDistinct(a)
		b := newBasis(a)

		// count must match 2^rank exactly.
		if uint64(len(b.d)) > kBits {
			t.Fatalf("rank overflow: a=%v rank=%d", a, len(b.d))
		}
		count := uint64(1) << uint(len(b.d))
		if count != uint64(len(d)) {
			t.Fatalf("count mismatch: a=%v basisCount=%d bruteCount=%d", a, count, len(d))
		}

		// every in-range k.
		for k := uint64(1); k <= count; k++ {
			got := b.kth(k)
			want := bruteKth(a, k)
			if got != want {
				t.Fatalf("mismatch: a=%v k=%d got=%d want=%d", a, k, got, want)
			}
		}

		// out-of-range k just past the count, and a couple further out -> -1.
		for _, k := range []uint64{count + 1, count + 2, count + 100} {
			if got := b.kth(k); got != -1 {
				t.Fatalf("oob: a=%v k=%d got=%d want=-1", a, k, got)
			}
		}
	}
}

// BenchmarkSolve measures the real-scale build+query envelope: n=1e5, q=1e5,
// 60-bit values. build is O(n*B), each query O(B).
func BenchmarkSolve(b *testing.B) {
	rng := rand.New(rand.NewPCG(42, 99))
	const n = 100000
	const q = 100000
	mask := (uint64(1) << kBits) - 1
	a := make([]uint64, n)
	for i := range a {
		a[i] = rng.Uint64() & mask
	}
	ks := make([]uint64, q)
	for i := range ks {
		ks[i] = 1 + rng.Uint64()%1000000000000000000 // 1..1e18
	}

	b.ResetTimer()
	var sink int64
	for i := 0; i < b.N; i++ {
		basis := newBasis(a)
		for _, k := range ks {
			sink ^= basis.kth(k)
		}
	}
	_ = sink
}
