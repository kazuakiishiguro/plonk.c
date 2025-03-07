#!/bin/bash

# Compile the unoptimized and optimized versions separately
echo "Compiling standard large benchmark..."
gcc -o src/large_benchmark_std.out src/large_benchmark.c -O0 -Wall

echo "Compiling optimized large benchmark..."
gcc -o src/large_benchmark_opt.out src/large_benchmark.c -O3 -march=native -flto -funroll-loops -Wall

# Run both benchmarks to compare
echo -e "\n======== STANDARD BUILD - LARGE CIRCUIT ========"
./src/large_benchmark_std.out

echo -e "\n======== OPTIMIZED BUILD - LARGE CIRCUIT ========"
./src/large_benchmark_opt.out

# Cleanup
rm src/large_benchmark_std.out src/large_benchmark_opt.out