// Minimal Chudnovsky implementation in Go (binary-splitting, arbitrary precision)
// Increase `n` (terms) or `digits` for more accuracy.
package main

import (
	"fmt"
	"math/big"
)

// binarySplit recursively returns (P,Q,R) for a≤k<b.
func binarySplit(a, b int64) (P, Q, R *big.Int) {
	if b == a+1 { // base case
		P = big.NewInt(0)
		Q = big.NewInt(0)
		R = big.NewInt(0)

		// P = −(6a−1)(2a−1)(6a−5)
		P.Mul(big.NewInt(6*a-1), big.NewInt(2*a-1))
		P.Mul(P, big.NewInt(6*a-5))
		P.Neg(P)

		// Q = 10939058860032000 · a³
		Q.Exp(big.NewInt(a), big.NewInt(3), nil)
		Q.Mul(Q, big.NewInt(10939058860032000))

		// R = P · (545140134 a + 13591409)
		R.Mul(P, big.NewInt(545140134*a+13591409))
		return
	}

	m := (a + b) / 2
	P1, Q1, R1 := binarySplit(a, m)
	P2, Q2, R2 := binarySplit(m, b)

	P = new(big.Int).Mul(P1, P2)
	Q = new(big.Int).Mul(Q1, Q2)
	R = new(big.Int).Add(
		new(big.Int).Mul(Q2, R1),
		new(big.Int).Mul(P1, R2),
	)
	return
}

// chudnovsky computes π using n terms at roughly `digits` decimal precision.
func chudnovsky(n int64, digits uint) *big.Float {
	precBits := digits * 4 // log₂10 ≈ 3.3 → 4× is safe
	_, Q, R := binarySplit(1, n)

	// coeff = 426880 · √10005
	sqrt10005 := new(big.Float).SetPrec(precBits).SetInt64(10005)
	sqrt10005.Sqrt(sqrt10005)
	coeff := new(big.Float).SetPrec(precBits).SetInt64(426880)
	coeff.Mul(coeff, sqrt10005)

	// numerator = coeff · Q
	num := new(big.Float).SetPrec(precBits).Mul(coeff, new(big.Float).SetInt(Q))

	// denominator = 13591409·Q + R
	den := new(big.Int).Mul(big.NewInt(13591409), Q)
	den.Add(den, R)

	return new(big.Float).SetPrec(precBits).Quo(num, new(big.Float).SetInt(den))
}

func main() {
	pi := chudnovsky(2, 80)       // 2 terms, ~80-digit precision
	fmt.Println(pi.Text('f', 30)) // print 30 decimals
}
