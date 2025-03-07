#!/bin/bash

# Compile and run the optimized benchmark
echo "=== Compiling Optimized Benchmark ==="

# Use advanced compiler flags for maximum optimization
CFLAGS="-O3 -march=native -mtune=native -flto -ffast-math -funroll-loops -Wall"

# Check if OpenMP is available
if [ -x "$(command -v brew)" ]; then
    # On macOS with Homebrew, we may need to specify paths for OpenMP
    if [ -d "$(brew --prefix libomp 2>/dev/null)" ]; then
        echo "Adding OpenMP support from Homebrew"
        OPENMP_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp"
        CFLAGS="$CFLAGS $OPENMP_FLAGS"
    fi
elif [ -x "$(command -v gcc)" ]; then
    # On Linux, GCC typically has OpenMP built-in
    CFLAGS="$CFLAGS -fopenmp"
fi

echo "Using compiler flags: $CFLAGS"

# Compile and run the benchmark
gcc $CFLAGS -o src/benchmark.out src/benchmark.c -lm

if [ $? -eq 0 ]; then
    echo "=== Running Optimized Benchmark ==="
    ./src/benchmark.out
else
    echo "Compilation failed."
fi

# Compile and run the optimization comparison benchmark
echo -e "\n=== Compiling Optimization Comparison Benchmark ==="
gcc $CFLAGS -o src/optimize_benchmark.out src/optimize_benchmark.c -lm

if [ $? -eq 0 ]; then
    echo "=== Running Optimization Comparison Benchmark ==="
    ./src/optimize_benchmark.out
else
    echo "Optimization comparison benchmark compilation failed."
fi

echo -e "\n=== Benchmark Complete ==="
echo "For detailed optimization information, see OPTIMIZATION_SUMMARY.md"