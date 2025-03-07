#!/bin/bash

# Compile the unoptimized and optimized versions separately
echo "Compiling standard benchmark..."
gcc -o src/benchmark_std.out src/benchmark.c -O0 -Wall

echo "Compiling optimized benchmark..."
gcc -o src/benchmark_opt.out src/benchmark.c -O3 -march=native -flto -funroll-loops -Wall

# Run both benchmarks to compare
echo -e "\n======== STANDARD BUILD ========"
./src/benchmark_std.out

echo -e "\n======== OPTIMIZED BUILD ========"
./src/benchmark_opt.out

# Cleanup
rm src/benchmark_std.out src/benchmark_opt.out