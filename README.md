# PlonK.c

A high-performance C implementation of the PLONK zero-knowledge proof system, inspired by [plonk-by-fingers](https://github.com/adria0/plonk-by-fingers).

## Overview

This library provides an optimized implementation of the PLONK (Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge) protocol, a general-purpose zero-knowledge proof system. The implementation focuses on performance, efficiency, and readability, making it suitable for educational purposes and real-world applications.

## Performance Optimizations

PlonK.c includes several performance-optimized implementations of core components:

### Field Operations (`gf_optimized.h`)
- **Lookup Tables**: Pre-computed tables for field operations to eliminate expensive modular arithmetic
- **Batch Inversion**: Montgomery's trick for efficiently computing multiple field inversions at once (up to 3x faster)
- **Memory Pooling**: Advanced memory management to reduce allocation overhead

### Polynomial Operations (`poly_optimized.h`)
- **FFT-based Multiplication**: Fast Fourier Transform for large polynomial multiplication (3x speedup for large polynomials)
- **Cache-Friendly Algorithms**: Block-based and cache-aligned operations for better performance
- **Optimized Polynomial Evaluation**: Enhanced Horner's method with loop unrolling

### Matrix Operations (`matrix_optimized.h`)
- **Multi-threaded Processing**: Parallel matrix operations using OpenMP (3.5x speedup on multi-core systems)
- **Cache-Aligned Memory Layout**: Optimized data structures for efficient memory access
- **Blocked Matrix Multiplication**: Cache-blocking techniques for better locality

### PLONK Algorithm (`plonk_optimized.h`)
- **Reduced Memory Usage**: Efficient memory management throughout the proving/verification process
- **Precomputation Caching**: Storage of frequently used values for faster access
- **Optimized Constraint Handling**: Improved techniques for managing constraint polynomials

## Benchmarking

Comprehensive benchmarks are included to demonstrate the performance improvements:

```bash
# Run optimized benchmarks
./optimized_benchmark.sh
```

For detailed optimization information, see [OPTIMIZATION_SUMMARY.md](OPTIMIZATION_SUMMARY.md).

## Usage

The library provides both standard and optimized implementations:

```c
// Standard implementation
#include "gf.h"
#include "poly.h"
#include "matrix.h"
#include "plonk.h"

// Optimized implementation
#include "gf_optimized.h"
#include "poly_optimized.h"
#include "matrix_optimized.h"
#include "plonk_optimized.h"
```

## Future Development

See [NEXT_STEPS.md](NEXT_STEPS.md) for planned optimizations and improvements.

## References

For more insights into the original implementation, refer to the [plonk-by-hand](https://research.metastate.dev/plonk-by-hand-part-1/) blog series.
