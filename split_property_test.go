package main

import (
	"math/big"
	"math/rand"
	"testing"
)

func eq(a, b *big.Int) bool {
	if a == nil || b == nil {
		return a == nil && b == nil
	}
	return a.Cmp(b) == 0
}

// TestParallelMatchesSerial checks that the parallel, FFT-using split produces
// bit-identical (P,Q,R) to the serial, stdlib-only reference across many
// ranges — including ranges large enough to exercise the FFT multiply path.
func TestParallelMatchesSerial(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	ranges := [][2]int64{
		{1, 2}, {1, 3}, {1, 10}, {2, 9}, {1, 100}, {50, 200}, {1, 1000},
	}
	for range 200 {
		a := int64(rng.Intn(2000) + 1)
		b := a + int64(rng.Intn(2000)+1)
		ranges = append(ranges, [2]int64{a, b})
	}
	// Large ranges to force operands past fftMinBits (~200k bits ≈ n>3800).
	if !testing.Short() {
		ranges = append(ranges, [2]int64{1, 8000}, [2]int64{1, 12000})
	}
	depth := parallelDepth()
	for _, r := range ranges {
		wP, wQ, wR := binarySplit(r[0], r[1])
		gP, gQ, gR := parallelSplit(r[0], r[1], depth, true)
		if !eq(wP, gP) || !eq(wQ, gQ) || !eq(wR, gR) {
			t.Fatalf("parallel != serial for [%d,%d)", r[0], r[1])
		}
	}
}

// TestRecurrence validates the binary-splitting combine formula at arbitrary
// (non-midpoint) split points.
func TestRecurrence(t *testing.T) {
	cases := [][3]int64{{1, 2, 10}, {1, 5, 10}, {1, 9, 10}, {3, 50, 200}, {1, 333, 1000}}
	for _, c := range cases {
		a, m, b := c[0], c[1], c[2]
		P, Q, R := binarySplit(a, b)
		P1, Q1, R1 := binarySplit(a, m)
		P2, Q2, R2 := binarySplit(m, b)
		cP := new(big.Int).Mul(P1, P2)
		cQ := new(big.Int).Mul(Q1, Q2)
		cR := new(big.Int).Add(new(big.Int).Mul(R1, Q2), new(big.Int).Mul(P1, R2))
		if !eq(P, cP) || !eq(Q, cQ) || !eq(R, cR) {
			t.Fatalf("recurrence broken at a=%d m=%d b=%d", a, m, b)
		}
	}
}

// TestSkipTopP confirms that skipping the top-level P (needP=false) leaves Q
// and R unchanged — only P is dropped.
func TestSkipTopP(t *testing.T) {
	depth := parallelDepth()
	for _, b := range []int64{10, 100, 1000, 8000} {
		_, Q1, R1 := parallelSplit(1, b, depth, true)
		Pn, Q0, R0 := parallelSplit(1, b, depth, false)
		if !eq(Q1, Q0) || !eq(R1, R0) {
			t.Fatalf("skip-top-P changed Q or R for n=%d", b)
		}
		// P is dropped (nil) only when the root recurses; below the serial
		// cutoff the leaf computes P and the caller simply ignores it.
		if b >= serialCutoff && Pn != nil {
			t.Fatalf("expected nil P when needP=false for n=%d", b)
		}
	}
}

// TestParallelDeterminism guards against goroutine-scheduling nondeterminism.
func TestParallelDeterminism(t *testing.T) {
	depth := parallelDepth()
	_, Q0, R0 := parallelSplit(1, 5000, depth, false)
	for i := range 20 {
		_, Q, R := parallelSplit(1, 5000, depth, false)
		if !eq(Q, Q0) || !eq(R, R0) {
			t.Fatalf("nondeterministic parallel split on run %d", i)
		}
	}
}

// TestMulMatchesStdlib exercises the FFT multiply path (including the negative
// operands that arise from P) against the standard library.
func TestMulMatchesStdlib(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	// Equal sizes, lopsided pairs, and pairs straddling the fftMinBits threshold.
	pairs := [][2]int{{1000, 1000}, {250000, 250000}, {500000, 500000},
		{500000, 1000}, {250000, 190000}, {210000, 190000}}
	for _, p := range pairs {
		x := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), uint(p[0])))
		y := new(big.Int).Rand(rng, new(big.Int).Lsh(big.NewInt(1), uint(p[1])))
		for _, sx := range []int64{1, -1} {
			for _, sy := range []int64{1, -1} {
				xs := new(big.Int).Mul(x, big.NewInt(sx))
				ys := new(big.Int).Mul(y, big.NewInt(sy))
				if !eq(mul(xs, ys), new(big.Int).Mul(xs, ys)) {
					t.Fatalf("mul mismatch at %d/%d bits, signs %d/%d", p[0], p[1], sx, sy)
				}
			}
		}
	}
}
