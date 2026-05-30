# chudnovsky

Computes π with the Chudnovsky algorithm — the whole expansion to arbitrary
precision, or a single digit deep in it. Binary splitting across all cores, and
FFT for the big multiplies, divides, and square root.

## Usage

```bash
go run . -digit 1000000     # the 1,000,000th digit
go run . -digit 100 -all    # π to 100 places
go run .                    # defaults to digit 10000
```

`-digit N` counts from the leading `3`, so `N=1` is `3` and `N=2` is the first
decimal. There's no base-10 spigot for π, so even a single deep digit means
computing the whole prefix — the digit mode just slices it.

## Speed

Total wall time to compute a digit, on a 10-core M-series:

| digit | time |
| --- | --- |
| 100,000 | 33 ms |
| 1,000,000 | 0.48 s |
| 10,000,000 | 9 s |

Go's `math/big` tops out at Karatsuba (`O(n^1.585)`) — no FFT — so the big
multiplies, divides, and √ go through
[bigfft](https://github.com/remyoudompheng/bigfft) (`O(n log n)`) above a size
threshold, with a stdlib fallback below it. The bigger the number, the more that
pays off.

## Notes

Pure Go, one dependency (bigfft, no cgo). MIT licensed.
