package main

import (
	"fmt"
	"testing"
)

func BenchmarkBinarySplit(b *testing.B) {
	for _, n := range []int64{100, 1000, 10000} {
		b.Run(fmt.Sprintf("n=%d", n), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				binarySplit(1, n)
			}
		})
	}
}

func BenchmarkParallelSplit(b *testing.B) {
	depth := parallelDepth()
	for _, n := range []int64{10000, 100000} {
		b.Run(fmt.Sprintf("n=%d", n), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				parallelSplit(1, n, depth, false)
			}
		})
	}
}

func BenchmarkExtractDigit(b *testing.B) {
	sizes := []int{1000, 10000, 100000}
	if !testing.Short() {
		sizes = append(sizes, 1000000)
	}
	for _, d := range sizes {
		b.Run(fmt.Sprintf("digit=%d", d), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				extractDigit(d)
			}
		})
	}
}
