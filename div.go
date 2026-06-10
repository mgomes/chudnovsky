package main

import "math/big"

var bigOne = big.NewInt(1)

// recip returns an approximation of ⌊2^s / v⌋ for v > 0, computed with Newton's
// method so the large multiplications go through the FFT path. Precision is
// doubled on the way down the recursion, so the cost is dominated by the final
// (full-size) step rather than O(log) full-size multiplies.
//
// The result is floor-biased: 0 ≤ 2^s/v − recip(v, s) < 2, independent of size.
// (Base case: exact floor, error in [0,1). Newton step: substituting the child
// error ε into r0·(2^(s+1) − v·r0) gives 2^(2s)/v − vε² exactly, and with
// sp = (s+vb)/2+1 the squared term is ≤ ε²/4 ulps, so the error recurrence
// E ← E²/4 + 1 has fixed point 2.)
func recip(v *big.Int, s uint) *big.Int {
	vb := uint(v.BitLen())
	// Base case: the result is small enough that stdlib division is cheaper than
	// another Newton level.
	if s <= vb+2*fftMinBits {
		num := new(big.Int).Lsh(bigOne, s)
		return num.Quo(num, v)
	}
	// Solve at roughly half precision, then refine with one Newton step:
	//   x ← x·(2^(s+1) − v·x) >> s   converges quadratically to 2^s/v.
	sp := (s+vb)/2 + 1
	r0 := new(big.Int).Lsh(recip(v, sp), s-sp) // ≈ 2^s/v at half precision
	vr := mul(v, r0)                           // ≈ 2^s
	t := new(big.Int).Lsh(bigOne, s+1)
	t.Sub(t, vr) // 2^(s+1) − v·r0
	r := mul(r0, t)
	return r.Rsh(r, s)
}

// divApprox returns ⌊u / v⌋ or one less, for u ≥ 0 and v > 0. With
// s = u.BitLen()+2 (the +2 is load-bearing), u·E/2^s < 1/2 for recip's
// E < 2, so the reciprocal multiply lands within one ulp below the exact
// floor. The π pipeline calls this directly: its guard digits absorb the
// slack, and skipping divFFT's verification multiply saves a full-size FFT.
// Small operands fall through to the standard library and are exact.
func divApprox(u, v *big.Int) *big.Int {
	if v.Sign() <= 0 {
		panic("divApprox: non-positive divisor")
	}
	if u.Sign() == 0 || u.Cmp(v) < 0 {
		return big.NewInt(0)
	}
	if u.BitLen() < 2*fftMinBits {
		return new(big.Int).Quo(u, v)
	}
	s := uint(u.BitLen() + 2)
	q := mul(u, recip(v, s)) // ≈ u·2^s/v
	return q.Rsh(q, s)
}

// divFFT returns exactly ⌊u / v⌋ for u ≥ 0 and v > 0: divApprox plus an exact
// correction. The remainder check costs a full-size multiply, so callers that
// can tolerate a one-ulp slack (the π pipeline) use divApprox directly.
func divFFT(u, v *big.Int) *big.Int {
	q := divApprox(u, v)
	if u.BitLen() < 2*fftMinBits || u.Cmp(v) < 0 {
		return q // those divApprox paths are exact
	}
	rem := new(big.Int).Sub(u, mul(q, v))
	for rem.Sign() < 0 {
		q.Sub(q, bigOne)
		rem.Add(rem, v)
	}
	for rem.Cmp(v) >= 0 {
		q.Add(q, bigOne)
		rem.Sub(rem, v)
	}
	return q
}
