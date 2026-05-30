// Parallel Chudnovsky implementation in Go (binary-splitting, arbitrary precision).
//
// Computes a specific decimal digit of π. Binary splitting runs across all CPU
// cores, and the large multiplications and the final division use FFT
// (github.com/remyoudompheng/bigfft) once the operands are big enough to beat
// the standard library's Karatsuba.
package main

import (
	"flag"
	"fmt"
	"math"
	"math/big"
	"runtime"
	"sync"
	"time"

	"github.com/remyoudompheng/bigfft"
)

// Per-term Chudnovsky constants, hoisted to avoid re-allocating them at every leaf.
var (
	cQBase     = big.NewInt(10939058860032000) // 640320³ / 24
	c13591409  = big.NewInt(13591409)
	c545140134 = int64(545140134)
)

const (
	// log2(10): bits required per decimal digit.
	log2of10 = 3.321928094887362
	// Chudnovsky yields log10(640320³/1728) ≈ 14.1816 decimal digits per term.
	digitsPerTerm = 14.181647462725477
	// Decimal guard digits computed beyond what is requested. The integer
	// pipeline is exact except for the floored √ and division, so the guard only
	// needs to exceed the longest run of 9s/0s near the target (≤ 6 below 10⁶,
	// well under 32 below ~5·10⁸).
	guardDigits = 32
	// Per-operand bit length at which bigfft.Mul beats stdlib Karatsuba on this
	// class of inputs (measured crossover ≈ 160k–200k bits; below it bigfft
	// switches to FFT too early and loses).
	fftMinBits = 200000
	// Context digits shown on each side of the target digit.
	ctxWindow = 5
	// Below this many terms a subtree is split serially — goroutine overhead
	// would otherwise dominate the tiny multiplies near the leaves.
	serialCutoff = 256
)

// mul returns x·y, using FFT for large operands and Karatsuba otherwise.
func mul(x, y *big.Int) *big.Int {
	if x.BitLen() >= fftMinBits && y.BitLen() >= fftMinBits {
		return bigfft.Mul(x, y)
	}
	return new(big.Int).Mul(x, y)
}

// pow10 returns 10^n by square-and-multiply, using the FFT path for the large
// products.
func pow10(n int) *big.Int {
	result := big.NewInt(1)
	base := big.NewInt(10)
	for n > 0 {
		if n&1 == 1 {
			result = mul(result, base)
		}
		if n >>= 1; n > 0 {
			base = mul(base, base)
		}
	}
	return result
}

// splitTerm is the binary-splitting base case for a single term k = a.
func splitTerm(a int64) (P, Q, R *big.Int) {
	// P = −(6a−1)(2a−1)(6a−5)
	P = new(big.Int).Mul(big.NewInt(6*a-1), big.NewInt(2*a-1))
	P.Mul(P, big.NewInt(6*a-5))
	P.Neg(P)

	// Q = (640320³/24)·a³
	Q = new(big.Int)
	if a <= 2_080_000 { // a³ < 2⁶³, so it fits in an int64
		Q.SetInt64(a * a * a)
	} else {
		Q.Exp(big.NewInt(a), big.NewInt(3), nil)
	}
	Q.Mul(Q, cQBase)

	// R = P·(545140134a + 13591409)
	R = new(big.Int).Mul(P, big.NewInt(c545140134*a+13591409))
	return
}

// binarySplit is the serial reference implementation over [a, b). It uses only
// the standard library multiply, so the tests can cross-check the parallel,
// FFT-using path against it.
func binarySplit(a, b int64) (P, Q, R *big.Int) {
	if b-a == 1 {
		return splitTerm(a)
	}
	m := (a + b) / 2
	P1, Q1, R1 := binarySplit(a, m)
	P2, Q2, R2 := binarySplit(m, b)
	P = new(big.Int).Mul(P1, P2)
	Q = new(big.Int).Mul(Q1, Q2)
	R = new(big.Int).Add(new(big.Int).Mul(R1, Q2), new(big.Int).Mul(P1, R2))
	return
}

// parallelSplit computes the binary split over [a, b) with the two halves and
// the combine multiplications run concurrently down to a serial cutoff.
//
// needP reports whether this node's P output is consumed by its parent. Only
// the left child's P feeds R = R1·Q2 + P1·R2, so the entire rightmost spine
// (starting at the root) can skip forming P = P1·P2 — the largest discarded
// multiply in the whole computation.
func parallelSplit(a, b int64, depth int, needP bool) (P, Q, R *big.Int) {
	if depth <= 0 || b-a < serialCutoff {
		P, Q, R = binarySplit(a, b)
		return
	}
	m := (a + b) / 2
	var P1, Q1, R1, P2, Q2, R2 *big.Int
	var wg sync.WaitGroup
	wg.Add(2)
	go func() { defer wg.Done(); P1, Q1, R1 = parallelSplit(a, m, depth-1, true) }()
	go func() { defer wg.Done(); P2, Q2, R2 = parallelSplit(m, b, depth-1, needP) }()
	wg.Wait()

	// Combine. Q, R1·Q2 and P1·R2 are independent products; run them
	// concurrently when the operands are large enough to be worth a goroutine.
	var qq, rq, pr, pp *big.Int
	if Q1.BitLen() >= fftMinBits {
		var cwg sync.WaitGroup
		run := func(dst **big.Int, x, y *big.Int) { defer cwg.Done(); *dst = mul(x, y) }
		cwg.Add(3)
		go run(&qq, Q1, Q2)
		go run(&rq, R1, Q2)
		go run(&pr, P1, R2)
		if needP {
			cwg.Add(1)
			go run(&pp, P1, P2)
		}
		cwg.Wait()
	} else {
		qq = mul(Q1, Q2)
		rq = mul(R1, Q2)
		pr = mul(P1, R2)
		if needP {
			pp = mul(P1, P2)
		}
	}
	Q = qq
	R = new(big.Int).Add(rq, pr)
	P = pp // nil when needP is false
	return
}

// terms returns the number of Chudnovsky terms needed for d correct decimals.
func terms(d int) int64 {
	return int64(math.Ceil(float64(d)/digitsPerTerm)) + 4
}

// parallelDepth picks the recursion depth that exposes ≈4 tasks per core.
func parallelDepth() int {
	d := 0
	for (1 << d) < 4*runtime.NumCPU() {
		d++
	}
	return d
}

// stageTimes records per-stage durations when piFloor is asked to profile.
type stageTimes struct{ split, sqrt, div time.Duration }

// piFloor returns ⌊π·10^d⌋ as a big.Int — its decimal string is "3" followed by
// the first d decimal digits of π. The whole pipeline is integer arithmetic
// (the only float is the one-off √10005), so the large multiply and the final
// division both go through the FFT path.
func piFloor(d int, st *stageTimes) *big.Int { return piFloorGuard(d, guardDigits, st) }

// piFloorGuard is piFloor with an explicit guard, so tests can confirm the
// result is stable as the guard grows (a too-small guard would diverge).
func piFloorGuard(d, guard int, st *stageTimes) *big.Int {
	total := d + guard

	t := time.Now()
	_, Q, R := parallelSplit(1, terms(total), parallelDepth(), false)
	if st != nil {
		st.split = time.Since(t)
	}

	// S = ⌊√10005 · 10^total⌋, via the FFT inverse-square-root.
	t = time.Now()
	S := sqrt10005Scaled(total)
	if st != nil {
		st.sqrt = time.Since(t)
	}

	// π·10^total = 426880·√10005·Q·10^total / (13591409·Q + R)
	//            = 426880·S·Q / (13591409·Q + R)
	t = time.Now()
	num := mul(new(big.Int).Mul(big.NewInt(426880), S), Q)
	den := new(big.Int).Add(new(big.Int).Mul(c13591409, Q), R)
	v := divFFT(num, den) // ⌊π·10^total⌋
	if st != nil {
		st.div = time.Since(t)
	}

	return v.Quo(v, pow10(guard)) // drop the guard digits → ⌊π·10^d⌋
}

// extractWindow returns the digit at digitPos (1-based; position 1 is the
// integer part '3') and an 11-character window centered on it.
func extractWindow(digitPos int, st *stageTimes) (digit int, window string) {
	v := piFloor(digitPos-1+ctxWindow, st) // materialize ctxWindow extra places
	width := 2*ctxWindow + 1
	win := new(big.Int).Mod(v, pow10(width))
	window = fmt.Sprintf("%0*d", width, win) // positions (digitPos-ctxWindow)..(digitPos+ctxWindow)
	digit = int(window[ctxWindow] - '0')
	return
}

// extractDigit returns just the digit at digitPos.
func extractDigit(digitPos int) int {
	d, _ := extractWindow(digitPos, nil)
	return d
}

// piDecimalString returns "3" followed by the first d decimal digits of π.
// Used by tests for full-string comparison; not on the hot path.
func piDecimalString(d int) string {
	if d == 0 {
		return "3"
	}
	return piFloor(d, nil).String()
}

func main() {
	digitPos := flag.Int("digit", 10000, "which digit of π to compute (digit 1 = the integer part '3')")
	verbose := flag.Bool("verbose", false, "print stage timings")
	flag.Parse()
	if *digitPos < 1 {
		*digitPos = 1
	}

	fmt.Printf("Using %d CPU cores\n", runtime.NumCPU())
	fmt.Printf("Calculating digit %d of π\n\n", *digitPos)

	var st stageTimes
	start := time.Now()
	digit, window := extractWindow(*digitPos, &st)
	elapsed := time.Since(start)

	fmt.Printf("Digit %d of π is: %d\n", *digitPos, digit)
	fmt.Printf("Total time: %v\n", elapsed)
	if *verbose {
		fmt.Printf("  split %v, sqrt %v, div %v\n", st.split, st.sqrt, st.div)
	}

	// Context window: trim positions before the integer part for small digitPos.
	before := window[:ctxWindow]
	after := window[ctxWindow+1:]
	if *digitPos <= ctxWindow {
		before = before[ctxWindow-(*digitPos-1):]
	}
	fmt.Printf("Context: ...%s[%d]%s...\n", before, digit, after)
}
