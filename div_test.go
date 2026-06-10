package main

import (
	"math/big"
	"math/rand"
	"testing"
)

// randBits returns a uniformly random integer of exactly `bits` bits.
func randBits(rng *rand.Rand, bits int) *big.Int {
	if bits <= 0 {
		return big.NewInt(0)
	}
	z := new(big.Int).Rand(rng, new(big.Int).Lsh(bigOne, uint(bits)))
	return z.SetBit(z, bits-1, 1) // ensure exactly `bits` bits
}

// TestDivFFT fuzzes the FFT divider against the standard library across sizes
// spanning the stdlib/FFT crossover, including lopsided and near-equal operands.
func TestDivFFT(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	pairs := [][2]int{
		{2048, 1024}, {500000, 250000}, {600000, 599000},
		{800000, 200000}, {1200000, 600000}, {2000000, 1000000},
	}
	for _, p := range pairs {
		u, v := randBits(rng, p[0]), randBits(rng, p[1])
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
	big1 := randBits(rng, 500000)
	if divFFT(big1, big1).Cmp(bigOne) != 0 {
		t.Fatal("u==v should be 1")
	}
	prod := mul(big1, randBits(rng, 500001))
	if divFFT(prod, big1).Cmp(new(big.Int).Quo(prod, big1)) != 0 {
		t.Fatal("exact multiple wrong")
	}
}

// TestRecipFloorBias asserts recip's error contract directly, as exact integer
// inequalities: 0 ≤ 2^s − v·r < 2v, i.e. the result is the floor of 2^s/v or
// one below. divApprox's one-ulp guarantee rests on this.
func TestRecipFloorBias(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	cases := []struct{ vbits, extra int }{
		{1024, 2048},     // stdlib base case
		{300000, 410000}, // one Newton level
		{600000, 900000}, // two Newton levels
		{599000, 1000},   // result far smaller than v
	}
	for _, c := range cases {
		v := randBits(rng, c.vbits)
		s := uint(c.vbits + c.extra)
		r := recip(v, s)
		diff := new(big.Int).Lsh(bigOne, s)
		diff.Sub(diff, mul(v, r)) // 2^s − v·r
		if diff.Sign() < 0 || diff.Cmp(new(big.Int).Lsh(v, 1)) >= 0 {
			t.Fatalf("recip out of contract for vbits=%d extra=%d: 2^s − v·r = %s ulps-ish",
				c.vbits, c.extra, new(big.Int).Quo(diff, v))
		}
	}
}

// TestDivApprox fuzzes the approximate divider the π pipeline uses: the result
// must be the exact floor or one below, never above.
func TestDivApprox(t *testing.T) {
	rng := rand.New(rand.NewSource(13))
	pairs := [][2]int{
		{2048, 1024}, {399000, 100000}, // stdlib path: exact
		{500000, 100000}, {900000, 450000}, // FFT reciprocal path
		{1200000, 600000}, {2000000, 1000000},
	}
	for _, p := range pairs {
		u, v := randBits(rng, p[0]), randBits(rng, p[1])
		q := divApprox(u, v)
		diff := new(big.Int).Quo(u, v)
		diff.Sub(diff, q) // ⌊u/v⌋ − q
		if diff.Sign() < 0 || diff.Cmp(bigOne) > 0 {
			t.Fatalf("divApprox off by %s for %d/%d bits (want 0 or 1 below the floor)", diff, p[0], p[1])
		}
	}
}
