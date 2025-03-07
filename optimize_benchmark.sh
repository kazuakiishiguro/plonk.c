#!/bin/bash

# Enable error reporting
set -e

# Compile the benchmark with advanced optimizations
echo "Compiling advanced optimization benchmark..."
gcc -o src/optimize_benchmark.out src/optimize_benchmark.c \
    -O3 -march=native -flto -ffast-math -funroll-loops -Wall \
    -I./src

# Run the benchmark
echo "Running advanced optimization benchmark..."
./src/optimize_benchmark.out

# Cleanup
rm src/optimize_benchmark.out

# Display summary
echo ""
echo "Performance summary of optimizations is available in OPTIMIZATION_SUMMARY.md"