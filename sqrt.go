package main

import (
	"math"
	"math/big"
)

// invSqrtBase is the precision (in bits) at which the Newton recursion bottoms
// out into a single stdlib integer √. Small enough that that √ is cheap, large
// enough to keep the recursion shallow.
const invSqrtBase = 1 << 16

// sqrt10005Scaled returns ⌊√10005 · 10^total⌋.
//
// It inverts the *small constant* 10005 with a Newton inverse-square-root, so
// the only large multiplications are the squaring y·y and the final ·10^total —
// both on ~total-digit operands on the FFT path (the c·y² term is a cheap
// scalar multiply). The whole thing is exact integer arithmetic — no big.Float,
// whose precision-tightening proved fragile at these sizes. Newton is
// floor-biased, so the result is within 1 of the exact floor; the caller's
// guard digits absorb that.
func sqrt10005Scaled(total int) *big.Int {
	p := uint(math.Ceil(float64(total)*log2of10)) + 64 // 2^p ≫ 10005·10^total
	r := invSqrtConst(10005, p)                        // ≈ ⌊2^p / √10005⌋
	s := new(big.Int).Mul(big.NewInt(10005), r)        // ≈ √10005 · 2^p
	s = mulPar(s, pow5(total))                         // · 5^total  (FFT)
	return s.Rsh(s, p-uint(total))                     // ⌊√10005 · 10^total⌋ — the ·2^total folds into the shift
}

// invSqrtConst returns ≈ ⌊2^p / √c⌋ for a small positive constant c, via Newton
// with precision doubling. The iteration y ← y·(3·2^(2p) − c·y²) >> (2p+1)
// converges quadratically. With t = 3·2^(2p) − c·y² written as 2^(2p+1) + δ,
// δ = 2^(2p) − c·y² (the Newton residual, ~1.5p bits), the step is computed —
// bit for bit identically — as y + ⌊y·δ/2^(2p+1)⌋, and both large products run
// on the unshifted half-precision ρ rather than on y = ρ·2^(p−pp), whose
// trailing zeros the FFT would otherwise pay for. Per-level multiply volume
// drops from ~5p to ~3p output bits.
func invSqrtConst(c int64, p uint) *big.Int {
	if p <= invSqrtBase {
		a := new(big.Int).Lsh(bigOne, 2*p)
		a.Quo(a, big.NewInt(c))
		return a.Sqrt(a)
	}
	pp := p/2 + 16
	rho := invSqrtConst(c, pp)       // ≈ 2^pp/√c, kept unshifted
	y2 := mulPar(rho, rho)           // ρ² (FFT squaring)
	cy2 := y2.Mul(y2, big.NewInt(c)) // scalar multiply
	delta := new(big.Int).Lsh(bigOne, 2*p)
	delta.Sub(delta, cy2.Lsh(cy2, 2*(p-pp))) // δ = 2^(2p) − c·y²  (y = ρ·2^(p−pp))
	yd := mulPar(rho, delta)                 // ρ·δ
	y := new(big.Int).Lsh(rho, p-pp)
	return y.Add(y, yd.Rsh(yd, p+pp+1)) // y + ⌊y·δ/2^(2p+1)⌋ = ⌊y·t/2^(2p+1)⌋
}
