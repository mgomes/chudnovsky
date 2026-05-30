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
	s = mul(s, pow10(total))                           // · 10^total  (FFT)
	return s.Rsh(s, p)                                 // ⌊√10005 · 10^total⌋
}

// invSqrtConst returns ≈ ⌊2^p / √c⌋ for a small positive constant c, via Newton
// with precision doubling. The iteration y ← y·(3·2^(2p) − c·y²) >> (2p+1)
// converges quadratically; c·y² is a scalar multiply, so y·y is the only large
// product.
func invSqrtConst(c int64, p uint) *big.Int {
	if p <= invSqrtBase {
		a := new(big.Int).Lsh(bigOne, 2*p)
		a.Quo(a, big.NewInt(c))
		return a.Sqrt(a)
	}
	pp := p/2 + 16
	y := new(big.Int).Lsh(invSqrtConst(c, pp), p-pp) // half-precision result, scaled up
	y2 := mul(y, y)                                  // FFT squaring (the only large product)
	cy2 := y2.Mul(y2, big.NewInt(c))                 // scalar multiply
	t := new(big.Int).Lsh(big.NewInt(3), 2*p)
	t.Sub(t, cy2) // 3·2^(2p) − c·y²
	y = mul(y, t)
	return y.Rsh(y, 2*p+1)
}
