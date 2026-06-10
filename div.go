package main

import "math/big"

var bigOne = big.NewInt(1)

// recipGuardBits is the relative-precision margin kept beyond a Newton level's
// result width when the divisor is truncated to that level's needs.
const recipGuardBits = 64

// recip returns an approximation of ⌊2^s / v⌋ for v > 0, computed with Newton's
// method so the large multiplications go through the FFT path. Precision is
// doubled on the way down the recursion, and the divisor is truncated to each
// level's working precision on the way, so the cost is dominated by the final
// (full-size) step rather than O(log) multiplies against the full-size v.
//
// Error contract: 2^s/v − recip(v, s) ∈ (−2^-50, 2) — within 2 below the exact
// value, or above it by an amount that can round the integer result one too
// high only when 2^s/v sits within 2^-50 of an integer. Below the bias bound:
// the Newton step is exact algebra, r0·(2^(s+1) − v·r0) = 2^(2s)/v − vε² for
// any child error ε, and with sp = (s+vb)/2+1 the vε² term is ≤ ε²/4 ulps, so
// the downward recurrence E ← E²/4 + 1 has fixed point 2. Above: truncating v
// to its top (s−vb)+64 bits (and lowering s by the same shift, which leaves
// the value unchanged) raises the target 2^(s−h)/⌊v/2^h⌋ relatively by less
// than 2^(2−64) per level, and the step never overshoots its own target.
func recip(v *big.Int, s uint) *big.Int {
	vb := uint(v.BitLen())
	// Bits of v below the level's result width (plus margin) cannot influence
	// the result beyond the contract — drop them so every multiply below, and
	// the base-case division, run on operands sized to the result.
	if s > vb && vb > (s-vb)+recipGuardBits {
		h := vb - (s - vb) - recipGuardBits
		v = new(big.Int).Rsh(v, h)
		s -= h
		vb -= h
	}
	// Base case: the result is small enough that stdlib division is cheaper than
	// another Newton level.
	if s <= vb+2*fftMinBits {
		num := new(big.Int).Lsh(bigOne, s)
		return num.Quo(num, v)
	}
	// Solve at roughly half precision, then refine with one Newton step,
	//   x ← x·(2^(s+1) − v·x) >> s,  quadratically convergent to 2^s/v,
	// computed — bit for bit identically — as r0 + ⌊ρ·d̃/2^(2sp−s)⌋ where
	// ρ is the unshifted child, r0 = ρ·2^(s−sp), and d̃ = 2^sp − v·ρ is the
	// child's residual (|d̃| < 2v). This keeps r0's s−sp trailing zero bits
	// and t's full width out of the FFTs.
	sp := (s+vb)/2 + 1
	rho := recip(v, sp) // ≈ 2^sp/v, kept unshifted
	d := new(big.Int).Lsh(bigOne, sp)
	d.Sub(d, mulPar(v, rho)) // d̃ = 2^sp − v·ρ
	rd := mulPar(rho, d)     // ρ·d̃
	r := new(big.Int).Lsh(rho, s-sp)
	return r.Add(r, rd.Rsh(rd, 2*sp-s)) // r0 + ⌊ρ·d̃/2^(2sp−s)⌋ = ⌊r0·t/2^s⌋
}

// divApprox returns ⌊u / v⌋ to within one ulp either side, for u ≥ 0 and
// v > 0. With s = u.BitLen()+2 (the +2 is load-bearing), u·E/2^s < 1/2 for
// recip's |E| < 2, so the reciprocal multiply lands on the exact floor, one
// below, or — only on a ~2^-50 fractional alignment — one above. The π
// pipeline calls this directly: its guard digits absorb the slack, and
// skipping divFFT's verification multiply saves a full-size FFT. Small
// operands fall through to the standard library and are exact.
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
	r := recip(v, s)
	// Bits of u below the quotient width (plus margin) shift the result by
	// well under an ulp; drop them so the final multiply is sized to the
	// quotient instead of to u.
	var k uint
	if vb := uint(v.BitLen()); vb > recipGuardBits {
		k = vb - recipGuardBits
	}
	q := mulPar(new(big.Int).Rsh(u, k), r) // ≈ (u/2^k)·(2^s/v)
	return q.Rsh(q, s-k)
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
