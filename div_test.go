package main

import (
	"math/big"
	"math/rand"
	"testing"
)

// TestDivFFT fuzzes the FFT divider against the standard library across sizes
// spanning the stdlib/FFT crossover, including lopsided and near-equal operands.
func TestDivFFT(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	randBits := func(bits int) *big.Int {
		if bits <= 0 {
			return big.NewInt(0)
		}
		z := new(big.Int).Rand(rng, new(big.Int).Lsh(bigOne, uint(bits)))
		return z.SetBit(z, bits-1, 1) // ensure exactly `bits` bits
	}
	pairs := [][2]int{
		{2048, 1024}, {500000, 250000}, {600000, 599000},
		{800000, 200000}, {1200000, 600000}, {2000000, 1000000},
	}
	for _, p := range pairs {
		u, v := randBits(p[0]), randBits(p[1])
		q := divFFT(u, v)
		if q.Cmp(new(big.Int).Quo(u, v)) != 0 {
			t.Fatalf("divFFT wrong for %d/%d bits", p[0], p[1])
		}
		// Postcondition both correction branches must establish: 0 ≤ u−q·v < v.
		// Asserting it directly exercises the invariant even though the
		// decrement branch is dormant under the current floor-biased recip.
		rem := new(big.Int).Sub(u, mul(q, v))
		if rem.Sign() < 0 || rem.Cmp(v) >= 0 {
			t.Fatalf("divFFT remainder out of range for %d/%d bits", p[0], p[1])
		}
	}
	// Edge cases: u<v, u==v, exact multiples.
	if divFFT(big.NewInt(3), big.NewInt(10)).Sign() != 0 {
		t.Fatal("u<v should be 0")
	}
	big1 := randBits(500000)
	if divFFT(big1, big1).Cmp(bigOne) != 0 {
		t.Fatal("u==v should be 1")
	}
	prod := mul(big1, randBits(500001))
	if divFFT(prod, big1).Cmp(new(big.Int).Quo(prod, big1)) != 0 {
		t.Fatal("exact multiple wrong")
	}
}
