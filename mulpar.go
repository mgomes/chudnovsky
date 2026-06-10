package main

import (
	"math/big"
	"math/bits"
	"sync"
)

// mulParMinBits gates the parallel multiply: below this per-operand size the
// goroutine and recombination overhead outweighs splitting the transform.
const mulParMinBits = 8 << 20

// mulPar returns x·y like mul, but splits large products into four quadrant
// multiplies run concurrently. bigfft's transform is single-threaded, so the
// serial stages (the √ Newton chain and the final division) otherwise occupy
// one core while the rest idle; FFT cost grows slightly faster than linearly
// in output size, so four half×half products finish in about half the wall
// time for ~2× the CPU. The quadrants slice the operands' limbs in place
// (read-only), and signs are reapplied at the end.
func mulPar(x, y *big.Int) *big.Int {
	if x.BitLen() < mulParMinBits || y.BitLen() < mulParMinBits {
		return mul(x, y)
	}
	xb, yb := x.Bits(), y.Bits()
	nx, ny := len(xb)/2, len(yb)/2
	var xl, xh, yl, yh big.Int
	xl.SetBits(xb[:nx])
	xh.SetBits(xb[nx:])
	yl.SetBits(yb[:ny])
	yh.SetBits(yb[ny:])

	var hh, hl, lh *big.Int
	var wg sync.WaitGroup
	wg.Add(3)
	go func() { defer wg.Done(); hh = mul(&xh, &yh) }()
	go func() { defer wg.Done(); hl = mul(&xh, &yl) }()
	go func() { defer wg.Done(); lh = mul(&xl, &yh) }()
	ll := mul(&xl, &yl)
	wg.Wait()

	const w = bits.UintSize // bits per big.Word
	z := new(big.Int).Lsh(hh, uint((nx+ny)*w))
	z.Add(z, hl.Lsh(hl, uint(nx*w)))
	z.Add(z, lh.Lsh(lh, uint(ny*w)))
	z.Add(z, ll)
	if (x.Sign() < 0) != (y.Sign() < 0) {
		z.Neg(z)
	}
	return z
}
