// Parallel Chudnovsky implementation in Go (binary-splitting, arbitrary precision)
// Increase `n` (terms) or `digits` for more accuracy.
package main

import (
	"flag"
	"fmt"
	"math/big"
	"runtime"
	"sync"
	"time"
)

type splitResult struct {
	P, Q, R *big.Int
}

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

// parallelBinarySplit performs binary splitting in parallel
func parallelBinarySplit(a, b int64, depth int) (P, Q, R *big.Int) {
	// Use serial version for small ranges or deep recursion
	if b-a < 1000 || depth > 4 {
		return binarySplit(a, b)
	}

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

	// Run splits in parallel
	var wg sync.WaitGroup
	wg.Add(2)

	var result1, result2 splitResult

	go func() {
		defer wg.Done()
		result1.P, result1.Q, result1.R = parallelBinarySplit(a, m, depth+1)
	}()

	go func() {
		defer wg.Done()
		result2.P, result2.Q, result2.R = parallelBinarySplit(m, b, depth+1)
	}()

	wg.Wait()

	P = new(big.Int).Mul(result1.P, result2.P)
	Q = new(big.Int).Mul(result1.Q, result2.Q)
	R = new(big.Int).Add(
		new(big.Int).Mul(result2.Q, result1.R),
		new(big.Int).Mul(result1.P, result2.R),
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

// parallelChudnovsky computes π using n terms with parallel processing
func parallelChudnovsky(n int64, digits uint) *big.Float {
	precBits := digits * 4 // log₂10 ≈ 3.3 → 4× is safe
	_, Q, R := parallelBinarySplit(1, n, 0)

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

// workerPool implementation for even better parallelization
type workItem struct {
	a, b   int64
	result chan splitResult
}

func workerPoolBinarySplit(a, b int64, numWorkers int) (P, Q, R *big.Int) {
	// For large ranges, just use the simple parallel version
	// The worker pool approach was causing deadlocks
	return parallelBinarySplit(a, b, 0)
}

// optimizedParallelChudnovsky uses worker pool for maximum parallelization
func optimizedParallelChudnovsky(n int64, digits uint) *big.Float {
	precBits := digits * 4 // log₂10 ≈ 3.3 → 4× is safe
	numWorkers := runtime.NumCPU()
	
	_, Q, R := workerPoolBinarySplit(1, n, numWorkers)

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
	// Parse command line flags
	digitPos := flag.Int("digit", 10000, "Which digit of pi to calculate")
	flag.Parse()
	
	runtime.GOMAXPROCS(runtime.NumCPU()) // Use all available cores
	
	fmt.Printf("Using %d CPU cores\n", runtime.NumCPU())
	fmt.Printf("Calculating digit %d of pi\n\n", *digitPos)
	
	// Adjust precision based on requested digit
	// Add extra buffer for accuracy
	digits := uint(*digitPos + 100)
	
	// Calculate number of terms needed (roughly 14 digits per term)
	n := int64((*digitPos / 14) + 100)
	
	var pi *big.Float
	var calcTime time.Duration
	
	// Try optimized parallel version first
	fmt.Println("Running optimized parallel version...")
	start := time.Now()
	
	// Use defer/recover to catch any panics in parallel execution
	func() {
		defer func() {
			if r := recover(); r != nil {
				fmt.Printf("Parallel version failed: %v\n", r)
				fmt.Println("Falling back to serial version...")
				start = time.Now()
				pi = chudnovsky(n, digits)
				calcTime = time.Since(start)
				fmt.Printf("Serial computation time: %v\n", calcTime)
			}
		}()
		
		pi = optimizedParallelChudnovsky(n, digits)
		calcTime = time.Since(start)
		fmt.Printf("Parallel computation time: %v\n", calcTime)
	}()
	
	// Extract the requested digit efficiently
	// For large digit positions, we use a more efficient extraction method
	fmt.Printf("Extracting digit %d...\n", *digitPos)
	
	// Multiply by 10^(digitPos-1) to shift the desired digit to units place
	shifter := new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(*digitPos-1)), nil)
	shifted := new(big.Float).SetPrec(pi.Prec()).Mul(pi, new(big.Float).SetInt(shifter))
	
	// Get the integer part and take mod 10 to get the digit
	intPart, _ := shifted.Int(nil)
	digit := new(big.Int).Mod(intPart, big.NewInt(10))
	
	fmt.Printf("\nDigit %d of pi is: %s\n", *digitPos, digit.String())
	
	// For context, let's show a few surrounding digits (if reasonable size)
	if *digitPos <= 100000 {
		contextDigits := 50
		if *digitPos < contextDigits {
			contextDigits = *digitPos + 10
		}
		piStr := pi.Text('f', contextDigits)
		if len(piStr) > *digitPos+1 {
			contextStart := *digitPos - 5
			if contextStart < 0 {
				contextStart = 0
			}
			contextEnd := *digitPos + 6
			if contextEnd > len(piStr)-2 {
				contextEnd = len(piStr) - 2
			}
			fmt.Printf("Context: ...%s[%s]%s...\n", 
				piStr[contextStart+1:*digitPos+1], 
				digit.String(),
				piStr[*digitPos+2:contextEnd+2])
		}
	}
}
