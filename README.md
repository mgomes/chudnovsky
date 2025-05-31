# Go Chudnovsky Algorithm

A high-performance parallel implementation of the Chudnovsky algorithm in Go for calculating π (pi) to arbitrary precision and extracting specific digits.

## Overview

This project provides a multi-threaded implementation of the [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm), one of the fastest methods for computing π. The implementation features:

- **Parallel binary splitting** for maximum performance
- **Efficient digit extraction** for specific positions
- **Multi-core utilization** using Go's goroutines
- **Arbitrary precision arithmetic** via Go's `math/big` package

## Features

- Calculate specific digits of π (e.g., the 1,000,000th digit)
- Parallel processing utilizing all CPU cores
- Optimized memory usage for large computations
- Fast execution (1.3 seconds for 1 million digits on modern hardware)

## Usage

Run the program with a command-line flag to specify which digit of π to calculate:

```bash
# Calculate the 1000th digit of π
go run main.go -digit 1000

# Calculate the 1,000,000th digit of π
go run main.go -digit 1000000

# Calculate the 10th digit of π (default if no flag provided)
go run main.go
```

### Example Output

```
Using 10 CPU cores
Calculating digit 1000000 of pi

Running optimized parallel version...
Parallel computation time: 1.345945708s
Extracting digit 1000000...

Digit 1000000 of pi is: 5
```

## The Chudnovsky Algorithm

The Chudnovsky algorithm is a fast method for calculating π based on Ramanujan's work. It converges much more rapidly than other methods, adding roughly 14 digits of π per term.

The formula used is:

```
1/π = 12 ∑ (-1)^k (6k)! (545140134k + 13591409) / ((3k)! (k!)^3 (640320^3)^k)
```

This implementation optimizes the computation using parallel binary splitting, which:
- Reduces asymptotic complexity from O(n²) to O(n log³ n)
- Utilizes multiple CPU cores for maximum performance
- Efficiently extracts specific digits without converting the entire number to string

## Performance

The parallel implementation provides significant speedup over serial computation:

- **1,000,000th digit**: ~1.3 seconds on modern multi-core systems
- **Scalable**: Performance scales with available CPU cores
- **Memory efficient**: Optimized digit extraction for large positions

## Architecture

The implementation includes several optimization levels:

1. **Serial binary splitting**: Basic optimized algorithm
2. **Parallel binary splitting**: Multi-threaded version with controlled recursion depth
3. **Optimized parallel version**: Production-ready implementation with efficient digit extraction

## Dependencies

This project only depends on Go's standard library:
- `flag` for command-line argument parsing
- `fmt` for output formatting
- `math/big` for arbitrary precision arithmetic
- `runtime` for CPU core detection
- `sync` for goroutine synchronization
- `time` for performance measurement

## License

[MIT License](LICENSE)
