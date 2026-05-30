# Go Chudnovsky Algorithm

A high-performance parallel implementation of the Chudnovsky algorithm in Go for
calculating π and extracting a specific decimal digit.

## Overview

This computes a chosen decimal digit of π via the
[Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm) with
binary splitting. To beat the standard library's Karatsuba ceiling at large
sizes, the heavy arithmetic uses FFT multiplication.

- **Parallel binary splitting** across all CPU cores
- **FFT multiplication** (Schönhage–Strassen, via
  [`bigfft`](https://github.com/remyoudompheng/bigfft)) for large operands, with
  an automatic fallback to `math/big` Karatsuba below the measured crossover
- **FFT-backed division and √** — a Newton reciprocal and a Newton
  inverse-square-root reduce both to FFT multiplies (the standard library's
  division and `Float.Sqrt` are Karatsuba-only)
- **Exact integer-domain pipeline** — the result is computed as `⌊π·10ⁿ⌋` with
  no floating point at all
- **Cheap digit extraction** via modular arithmetic (no full base conversion)

## Usage

`-digit N` selects the digit by 1-based position, where **position 1 is the
integer part `3`** and position N (N ≥ 2) is the (N−1)th decimal digit. So the
famous millionth *decimal* digit (`1`) is `-digit 1000001`.

```bash
go run . -digit 1000          # the 1000th position
go run . -digit 1000000       # the 1,000,000th position
go run . -digit 1000000 -verbose   # also print per-stage timings
go run .                      # default: digit 10000
```

### Example output

```
Using 10 CPU cores
Calculating digit 1000000 of π

Digit 1000000 of π is: 5
Total time: 545ms
Context: ...94581[5]13092...
```

## The algorithm

Chudnovsky converges at ≈14.18 digits per term:

```
1/π = 12 ∑ (-1)^k (6k)! (545140134k + 13591409) / ((3k)! (k!)³ (640320³)^k)
```

Binary splitting evaluates the truncated series as a single rational `Q`, `R`,
turning the sum into a balanced tree of big-integer multiplications — `O(M(n)
log n)` instead of `O(n²)`. The largest of those multiplications (and the final
division) dominate at scale, which is where FFT pays off.

## Performance

Measured on a 10-core Apple Silicon machine, total wall time (best of N):

| digit position | before | after | speedup |
| --- | --- | --- | --- |
| 10,000 | 40 ms | 5 ms | 7.6× |
| 100,000 | 3.5 s | 33 ms | 105× |
| 1,000,000 | 1.32 s | 0.48 s | 2.8× |
| 10,000,000 | 52.6 s | 9.0 s | 5.9× |

The 10k–100k jumps come largely from fixing an extraction bug; the 1M–10M gains
are the FFT + parallel arithmetic, and grow with size. Note that Go's `math/big`
has no FFT multiply, so the standard-library ceiling is Karatsuba (`O(n^1.585)`);
`bigfft` brings the hot multiply/divide/√ down toward `O(n log n)`.

## Tests

```bash
go test -short -race ./...   # fast: unit + property tests, race detector
go test -race ./...          # full: includes the 1,000,000-digit regression lock
go test -bench=. ./...       # benchmarks (binary split, parallel split, extraction)
```

The suite locks correctness against a reference value of π (1000 decimals),
checks the parallel/FFT path bit-for-bit against the serial Karatsuba reference,
fuzzes the FFT divider against `math/big`, and pins documented digits as
regression locks.

## Dependencies

- Go's standard library (`math/big`, `runtime`, `sync`, …)
- [`github.com/remyoudompheng/bigfft`](https://github.com/remyoudompheng/bigfft)
  — pure-Go FFT big-integer multiplication (no cgo)

## License

[MIT License](LICENSE)
