# chudnovsky

A high-performance parallel implementation of the Chudnovsky algorithm in Go.
Computes π to arbitrary precision and extracts a specific decimal digit.

## Overview

This computes π to arbitrary precision — the full expansion (`-all`) or a single
chosen digit deep in it — via the
[Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm) with
binary splitting. To beat the standard library's Karatsuba ceiling at large
sizes, the heavy arithmetic uses FFT multiplication. (There is no base-10
spigot, so a single deep digit still requires computing the whole prefix — the
digit mode just slices the same arbitrary-precision value.)

- **Parallel binary splitting** across all CPU cores
- **FFT multiplication** (Schönhage–Strassen, via
  [`bigfft`](https://github.com/remyoudompheng/bigfft)) for large operands, with
  an automatic fallback to `math/big` Karatsuba below the measured crossover
- **Truncated-Newton division and √** — a Newton reciprocal and a Newton
  inverse-square-root reduce both to FFT multiplies (the standard library's
  division and `Float.Sqrt` are Karatsuba-only), every operand is truncated to
  the precision its step actually needs, and the Newton steps multiply only the
  half-width residuals
- **Integer-domain pipeline** — the result is computed as `⌊π·10ⁿ⌋` with no
  floating point at all; the few bounded approximations (floored √, reciprocal
  slack, operand truncation) stay within ulps of `10ⁿ⁺ᵍ` and are absorbed by
  the guard digits
- **Quadrant-parallel large multiplies** keep the serial √/division stages on
  all cores (bigfft's transform is single-threaded)
- **Cheap digit extraction** via modular arithmetic (no full base conversion)

## Usage

`-digit N` selects the digit by 1-based position, where **position 1 is the
integer part `3`** and position N (N ≥ 2) is the (N−1)th decimal digit. So the
famous millionth *decimal* digit (`1`) is `-digit 1000001`.

```bash
go run . -digit 1000          # the 1000th position
go run . -digit 1000000       # the 1,000,000th position
go run . -digit 100 -all      # print π to 100 places
go run . -digit 1000000 -verbose   # also print per-stage timings
go run .                      # default: digit 10000
```

`-all` prints the full expansion to `-digit` places instead of just the digit at
that position.

### Example output

```
Using 10 CPU cores
Calculating digit 1000000 of π

Digit 1000000 of π is: 5
Total time: 194ms
Context: ...94581[5]13092...
```

## The algorithm

Chudnovsky converges at ≈14.18 digits per term:

```math
\frac{1}{\pi}
=
\frac{12}{640320^{3/2}}
\sum_{k=0}^{\infty}
\frac{(-1)^k (6k)! (13591409 + 545140134k)}
{(3k)! (k!)^3 640320^{3k}}.
```

Equivalently, the implementation truncates the denominator series:

```math
S_n =
\sum_{k=0}^{n-1}
\frac{(-1)^k (6k)! (13591409 + 545140134k)}
{(3k)! (k!)^3 640320^{3k}},
\qquad
\pi \approx \frac{426880\sqrt{10005}}{S_n}.
```

Binary splitting evaluates the truncated series as a single rational `Q`, `R`,
turning the sum into a balanced tree of big-integer multiplications — `O(M(n)
log n)` instead of `O(n²)`. The largest of those multiplications (and the final
division) dominate at scale, which is where FFT pays off.

Exact splitting leaves `Q` and `R` ≈2.3× larger than the final quotient needs,
and the quotient `426880·S·Q / (13591409·Q + R)` is invariant when both scale
together — so both are truncated to the target precision (plus guard) before
the division, and the Newton reciprocal behind that division re-truncates its
divisor at every precision-doubling level. The √ runs concurrently with the
splitting, which hides it entirely.

## Performance

Measured on a 10-core Apple Silicon machine, total wall time. The round 1 and
round 2 columns were re-measured back to back on the same machine and Go
toolchain (best of 3, interleaved); the serial baseline is the original
pre-optimization implementation.

| digit position | serial baseline | round 1: FFT + parallel | round 2: truncated Newton | total |
| --- | --- | --- | --- | --- |
| 10,000 | 40 ms | 2.6 ms | 2.0 ms | 20× |
| 100,000 | 3.5 s | 32 ms | 17 ms | 206× |
| 1,000,000 | 1.32 s | 0.51 s | 0.19 s | 6.9× |
| 10,000,000 | 52.6 s | 9.5 s | 2.2 s | 24× |

Round 1 (the 10k–100k jumps come largely from fixing an extraction bug) routed
the heavy arithmetic through FFT and parallelized the splitting: Go's
`math/big` has no FFT multiply, so the standard-library ceiling is Karatsuba
(`O(n^1.585)`); `bigfft` brings the hot multiply/divide/√ down toward
`O(n log n)`. Round 2 stopped doing work the result can't see: at 10M digits
the division stage was 69% of the runtime, dividing numbers ≈3.3× wider than
the quotient with a reciprocal computed against the full-width divisor at
every Newton level. Truncating to needed precision, multiplying only Newton
residuals, overlapping the √ with the splitting, and running the remaining
serial multiplies in quadrants took 10M digits from 9.5s to 2.2s (div stage:
6.6s → 0.7s).

## Tests

```bash
go test -short -race ./...   # fast: unit + property tests, race detector
go test -race ./...          # full: includes the 1,000,000-digit regression lock
go test -bench=. ./...       # benchmarks (binary split, parallel split, extraction)
```

The suite locks correctness against a reference value of π (1000 decimals),
checks the parallel/FFT path bit-for-bit against the serial Karatsuba reference,
fuzzes the FFT divider and inverse-square-root against `math/big`, and pins
documented digits as regression locks.

## Dependencies

- Go's standard library (`math/big`, `runtime`, `sync`, …)
- [`github.com/remyoudompheng/bigfft`](https://github.com/remyoudompheng/bigfft)
  — pure-Go FFT big-integer multiplication (no cgo)

## License

[MIT License](LICENSE)
