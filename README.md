# chudnovsky

Computes π with the Chudnovsky algorithm — the whole expansion to arbitrary
precision, or a single digit deep in it. Binary splitting runs across all cores;
FFT, truncated Newton division, and a Newton inverse square root handle the big
integer work.

## Usage

```bash
go run . -digit 1000              # the 1000th position
go run . -digit 1000000           # the 1,000,000th position
go run . -digit 100 -all          # π to 100 decimal places
go run . -digit 1000000 -verbose  # also print per-stage timings
go run .                          # defaults to digit 10000
```

`-digit N` counts from the leading `3`: `N=1` is `3`, and `N=2` is the first
decimal. The famous millionth decimal digit is therefore `-digit 1000001`.
There's no base-10 spigot for π, so even a single deep digit means computing the
whole prefix; digit mode just slices it.

## Speed

Total wall time to compute a digit, on a 10-core M-series:

| digit position | time |
| --- | --- |
| 10,000 | 2.0 ms |
| 100,000 | 17 ms |
| 1,000,000 | 0.19 s |
| 10,000,000 | 2.2 s |

Go's `math/big` tops out at Karatsuba (`O(n^1.585)`) — no FFT — so the big
multiplies, divides, and √ go through
[bigfft](https://github.com/remyoudompheng/bigfft) (`O(n log n)`) above a size
threshold, with a stdlib fallback below it. The bigger the number, the more that
pays off.

## Notes

```
1/π = 12 ∑ (-1)^k (6k)! (545140134k + 13591409) / ((3k)! (k!)³ (640320³)^k)
```

Binary splitting evaluates the truncated series as a single rational `Q`, `R`,
turning the sum into a balanced tree of big-integer multiplications. Exact
splitting leaves `Q` and `R` wider than the final quotient needs, so both are
truncated to the target precision plus guard digits before division. The
division uses a truncated Newton reciprocal, √10005 uses a Newton inverse square
root, and the remaining serial large multiplies split into quadrants.

The result is computed as `⌊π·10ⁿ⌋` with integer arithmetic. The bounded
approximations from the floored √, reciprocal slack, and operand truncation are
absorbed by guard digits.

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

Pure Go, one dependency ([bigfft](https://github.com/remyoudompheng/bigfft), no
cgo). MIT licensed.
