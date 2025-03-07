#!/bin/bash

# Compile benchmark with high optimization
echo "Compiling optimized benchmark..."
gcc -o src/benchmark.out src/benchmark.c -O3 -march=native -flto -funroll-loops -Wall

# Run benchmark
echo "Running optimized benchmark..."
./src/benchmark.out

# Cleanup
rm src/benchmark.out