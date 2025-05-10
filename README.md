# Go Chudnovsky Algorithm

A minimal and efficient implementation of the Chudnovsky algorithm in Go for calculating π (pi) to arbitrary precision.

## Overview

This project provides a high-performance implementation of the [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm), one of the fastest methods for computing π. The implementation uses:

- Binary splitting technique for improved performance
- Arbitrary precision arithmetic via Go's `math/big` package

## The Chudnovsky Algorithm

The Chudnovsky algorithm is a fast method for calculating π based on Ramanujan's work. It converges much more rapidly than other methods, adding roughly 14 digits of π per term.

The formula used is:

```
1/π = 12 ∑ (-1)^k (6k)! (545140134k + 13591409) / ((3k)! (k!)^3 (640320^3)^k)
```

This implementation optimizes the computation using binary splitting, which reduces the asymptotic complexity from O(n²) to O(n log³ n).

## Usage

```go
package main

import (
    "fmt"
    "github.com/mgomes/go-chudnovsky"
)

func main() {
    // Calculate π with 2 terms (~80 decimal digits precision)
    pi := chudnovsky(2, 80)

    // Print the first 30 decimal digits
    fmt.Println(pi.Text('f', 30))
}
```

### Adjusting Precision

You can modify two main parameters to adjust the precision:

1. `n` (terms): The number of terms in the series to compute (each term adds ~14 digits)
2. `digits`: The desired decimal precision

Examples:

```go
// ~120 digits of precision (3 terms)
pi := chudnovsky(3, 120)

// ~1000 digits of precision (100 terms)
pi := chudnovsky(100, 1000)
```

## Performance

The implementation uses binary splitting, which significantly outperforms naive summation for large numbers of terms. This makes it practical to compute thousands or even millions of digits of π.

## Dependencies

This project only depends on Go's standard library:
- `fmt` for output formatting
- `math/big` for arbitrary precision arithmetic

## License

[MIT License](LICENSE)
