package main

import (
	"math/big"
	"testing"
)

// TestSqrt10005Scaled checks the FFT inverse-square-root against the standard
// library's exact integer √ of 10005·10^(2·total), across sizes that exercise
// the base case and several Newton levels. Newton is floor-biased, so the
// result must land in [exact−1, exact]; the pi pipeline's guard digits absorb
// that sub-ulp slack.
func TestSqrt10005Scaled(t *testing.T) {
	totals := []int{1000, 20000, 70000, 130000, 300000}
	if !testing.Short() {
		totals = append(totals, 1000000)
	}
	for _, total := range totals {
		got := sqrt10005Scaled(total)
		exact := new(big.Int).Sqrt(new(big.Int).Mul(big.NewInt(10005), pow10(2*total)))
		diff := new(big.Int).Sub(exact, got)
		if diff.Sign() < 0 || diff.Cmp(big.NewInt(1)) > 0 {
			t.Fatalf("sqrt10005Scaled(%d) off by %s (must be 0 or 1 below exact)", total, diff)
		}
	}
}
