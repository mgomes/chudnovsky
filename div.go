package main

import "math/big"

var bigOne = big.NewInt(1)

// recip returns an approximation of ⌊2^s / v⌋ for v > 0, computed with Newton's
// method so the large multiplications go through the FFT path. Precision is
// doubled on the way down the recursion, so the cost is dominated by the final
// (full-size) step rather than O(log) full-size multiplies.
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
	vr := mul(v, r0)                            // ≈ 2^s
	t := new(big.Int).Lsh(bigOne, s+1)
	t.Sub(t, vr) // 2^(s+1) − v·r0
	r := mul(r0, t)
	return r.Rsh(r, s)
}

// divFFT returns ⌊u / v⌋ for u ≥ 0 and v > 0, using FFT multiplication for the
// large products. It falls back to the standard library when the operands are
// too small for FFT to pay off.
func divFFT(u, v *big.Int) *big.Int {
	if v.Sign() <= 0 {
		panic("divFFT: non-positive divisor")
	}
	if u.Sign() == 0 || u.Cmp(v) < 0 {
		return big.NewInt(0)
	}
	if u.BitLen() < 2*fftMinBits {
		return new(big.Int).Quo(u, v)
	}
	s := uint(u.BitLen() + 2)
	q := mul(u, recip(v, s)) // ≈ u·2^s/v
	q.Rsh(q, s)              // ≈ u/v, accurate to within a few ulps
	// Exact correction: Newton + floor can leave q off by a small amount.
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
