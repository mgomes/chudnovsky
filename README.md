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

10-core M-series, total wall time:

| digit | before | after | |
| --- | --- | --- | --- |
| 100,000 | 3.5 s | 33 ms | 105× |
| 1,000,000 | 1.3 s | 0.48 s | 2.8× |
| 10,000,000 | 53 s | 9 s | 5.9× |

Go's `math/big` has no FFT, so it tops out at Karatsuba (`O(n^1.585)`). Dropping
in [bigfft](https://github.com/remyoudompheng/bigfft) for the hot
multiply/divide/√ pulls it toward `O(n log n)`, and the win grows with size. (The
100k jump is mostly a separate fix — it used to spend seconds formatting a few
context digits.)

## Notes

Pure Go, one dependency (bigfft, no cgo). MIT licensed.
