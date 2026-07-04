package main

import (
	"bufio"
	"math/bits"
	"os"
	"strconv"
)

// k-th smallest distinct subset-XOR.
//
// the set of every XOR you can build from a subset of the array is a linear
// subspace over GF(2). the array elements span it. the empty subset gives 0, so
// 0 is always in -- that's the zero vector, and every subspace has it.
//
// the honest way is to enumerate 2^n subsets. dead at n = 1e5. so stop thinking
// about subsets and think about the subspace they span.
//
// build a reduced basis. r = rank = number of basis vectors. a subspace of
// dimension r over GF(2) has exactly 2^r points -- each of the r basis vectors is
// either in your combination or not, and distinct choices give distinct vectors
// because the basis is independent. so the count of distinct achievable XORs is
// 2^r. no duplicates, no enumeration.
//
// the ordering is the trick. reduce the basis to row-echelon-and-back: each basis
// vector owns a unique highest bit, and that owned bit is cleared from every other
// basis vector. sort the r vectors ascending -- ascending by value is the same as
// ascending by owned bit, because each is the smallest vector carrying its top
// bit. now the j-th basis vector (0-indexed, ascending) behaves like the j-th
// binary digit of the answer's rank. why: its owned bit sits in no other basis
// vector, so flipping it in or out moves the value across a clean bit boundary
// that no lower basis vector can reach. the rank-(k-1) point is the XOR of the
// basis vectors whose position bit is set in (k-1).
//
// lower bound. you must read all n inputs -- Omega(n). you must answer all q
// queries -- Omega(q). basis build touches each element across B bits: O(n*B).
// each query scans r <= B bits: O(q*B). B = 60 here. so O((n+q)*B) total, which
// is Omega(n+q) up to the constant B -- optimal to the width of the numbers.
const kBits = 60 // a_i < 2^60, so bit positions 0..59, rank at most 60.

// xorBasis holds the reduced basis vectors, ascending by owned highest bit.
// d[j] is the j-th smallest -- exactly the j-th binary digit of a query rank.
type xorBasis struct {
	d []uint64
}

// newBasis builds the reduced GF(2) basis of a. O(len(a) * kBits).
func newBasis(a []uint64) xorBasis {
	// p[b] holds the basis vector that owns highest bit b, or 0 if bit b is free.
	var p [kBits]uint64
	for _, x := range a {
		// fold x down through owned pivots until it finds a free bit or dies.
		for x != 0 {
			b := bits.Len64(x) - 1 // highest set bit of x
			if p[b] == 0 {
				p[b] = x // x owns bit b now
				break
			}
			x ^= p[b] // clear bit b, keep folding
		}
		// x == 0 means it was already spanned -- a dependent vector, dropped.
	}

	// reduce to full row-echelon: strip every owned bit out of the higher vectors.
	// after this each p[i] carries pivot i and no other pivot -- so ascending index
	// is ascending value, and each digit of a rank flips one clean bit boundary.
	for i := 0; i < kBits; i++ {
		if p[i] == 0 {
			continue
		}
		for j := 0; j < i; j++ {
			if p[j] != 0 && p[i]&(uint64(1)<<j) != 0 {
				p[i] ^= p[j]
			}
		}
	}

	// collect the live vectors low-to-high -- that's ascending order, ready to index.
	b := xorBasis{d: make([]uint64, 0, kBits)}
	for i := 0; i < kBits; i++ {
		if p[i] != 0 {
			b.d = append(b.d, p[i])
		}
	}
	return b
}

// kth returns the k-th smallest distinct achievable XOR (k is 1-indexed), or -1
// when k runs past the 2^r values that exist. O(rank).
func (b xorBasis) kth(k uint64) int64 {
	r := len(b.d)
	// count = 2^r. r <= kBits = 60, so 1<<r fits a uint64 with room -- no overflow,
	// no clever compare needed. guard the shift anyway in case kBits ever grows.
	if r < 64 {
		if k > (uint64(1) << uint(r)) {
			return -1 // asked past the last achievable value
		}
	}
	idx := k - 1 // 0-indexed rank
	var x uint64
	// idx < 2^r, so only bits 0..r-1 can be set -- each picks its basis vector.
	for j := 0; j < r; j++ {
		if idx&(uint64(1)<<uint(j)) != 0 {
			x ^= b.d[j]
		}
	}
	return int64(x)
}

func main() {
	in := bufio.NewReaderSize(os.Stdin, 1<<20)
	out := bufio.NewWriterSize(os.Stdout, 1<<20)
	defer out.Flush()

	readUint := func() uint64 {
		var v uint64
		c, err := in.ReadByte()
		// skip anything that isn't a digit -- whitespace, newlines, EOF slack.
		for err == nil && (c < '0' || c > '9') {
			c, err = in.ReadByte()
		}
		for err == nil && c >= '0' && c <= '9' {
			v = v*10 + uint64(c-'0')
			c, err = in.ReadByte()
		}
		return v
	}

	n := int(readUint())
	a := make([]uint64, n)
	for i := range a {
		a[i] = readUint()
	}
	basis := newBasis(a)

	q := int(readUint())
	buf := make([]byte, 0, 24)
	for ; q > 0; q-- {
		k := readUint()
		buf = strconv.AppendInt(buf[:0], basis.kth(k), 10)
		buf = append(buf, '\n')
		out.Write(buf)
	}
}
