#!/bin/bash

# Compile benchmark
echo "Compiling benchmark..."
gcc -o src/benchmark.out src/benchmark.c -O3 -Wall

# Run benchmark
echo "Running benchmark..."
./src/benchmark.out

# Cleanup
rm src/benchmark.out