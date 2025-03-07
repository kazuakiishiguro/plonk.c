#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "gf.h"

// Timer function with microsecond precision
double get_time_us() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000.0) + tv.tv_usec;
}

// Number of iterations for benchmarks
#define BENCHMARK_ITERATIONS 100000

// Precomputed tables for optimized operations
static uint8_t GF_ADD_TABLE[MODULO_GF][MODULO_GF];
static uint8_t GF_MUL_TABLE[MODULO_GF][MODULO_GF];
static uint8_t GF_INV_TABLE[MODULO_GF];
static bool TABLES_INITIALIZED = false;

// Initialize tables for optimized operations
void init_tables() {
    if (TABLES_INITIALIZED) return;
    
    // Addition table
    for (int i = 0; i < MODULO_GF; i++) {
        for (int j = 0; j < MODULO_GF; j++) {
            GF_ADD_TABLE[i][j] = (i + j) % MODULO_GF;
        }
    }
    
    // Multiplication table
    for (int i = 0; i < MODULO_GF; i++) {
        for (int j = 0; j < MODULO_GF; j++) {
            GF_MUL_TABLE[i][j] = (i * j) % MODULO_GF;
        }
    }
    
    // Inverse table using Fermat's Little Theorem
    GF_INV_TABLE[0] = 0; // Special case: 0 has no inverse
    for (int i = 1; i < MODULO_GF; i++) {
        int result = 1;
        int base = i;
        int exp = MODULO_GF - 2;
        
        while (exp > 0) {
            if (exp & 1) result = (result * base) % MODULO_GF;
            base = (base * base) % MODULO_GF;
            exp >>= 1;
        }
        
        GF_INV_TABLE[i] = result;
    }
    
    TABLES_INITIALIZED = true;
}

// Optimized field operations with lookup tables
GF opt_gf_add(GF a, GF b) {
    GF result = { GF_ADD_TABLE[a.value][b.value] };
    return result;
}

GF opt_gf_mul(GF a, GF b) {
    GF result = { GF_MUL_TABLE[a.value][b.value] };
    return result;
}

GF opt_gf_inv(GF a) {
    GF result = { GF_INV_TABLE[a.value] };
    return result;
}

// Benchmark field operations
void benchmark_field_operations() {
    printf("Benchmarking field operations (%d iterations)...\n", BENCHMARK_ITERATIONS);
    
    // Create random field elements
    GF* a = (GF*)malloc(BENCHMARK_ITERATIONS * sizeof(GF));
    GF* b = (GF*)malloc(BENCHMARK_ITERATIONS * sizeof(GF));
    
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        a[i] = gf_new(rand() % MODULO_GF);
        b[i] = gf_new(rand() % MODULO_GF);
    }
    
    // Original implementation
    {
        // Test addition
        double start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
            GF result = gf_add(a[i], b[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        double end = get_time_us();
        printf("  Original GF add: %.2f us\n", (end - start) / BENCHMARK_ITERATIONS);
        
        // Test multiplication
        start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
            GF result = gf_mul(a[i], b[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        end = get_time_us();
        printf("  Original GF mul: %.2f us\n", (end - start) / BENCHMARK_ITERATIONS);
        
        // Test inversion (expensive)
        start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS / 100; i++) {
            GF result = gf_inv(a[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        end = get_time_us();
        printf("  Original GF inv: %.2f us\n", (end - start) / (BENCHMARK_ITERATIONS / 100));
    }
    
    printf("\n");
    
    // Optimized implementation
    {
        // Initialize tables
        init_tables();
        
        // Test addition (table-based)
        double start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
            GF result = opt_gf_add(a[i], b[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        double end = get_time_us();
        printf("  Optimized GF add: %.2f us\n", (end - start) / BENCHMARK_ITERATIONS);
        
        // Test multiplication (table-based)
        start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
            GF result = opt_gf_mul(a[i], b[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        end = get_time_us();
        printf("  Optimized GF mul: %.2f us\n", (end - start) / BENCHMARK_ITERATIONS);
        
        // Test inversion (table-based)
        start = get_time_us();
        for (int i = 0; i < BENCHMARK_ITERATIONS / 100; i++) {
            GF result = opt_gf_inv(a[i]);
            if (result.value == 0) printf("."); // Prevent optimization
        }
        end = get_time_us();
        printf("  Optimized GF inv: %.2f us\n", (end - start) / (BENCHMARK_ITERATIONS / 100));
    }
    
    free(a);
    free(b);
}

int main() {
    // Seed random number generator
    srand(time(NULL));
    
    printf("=== Advanced Optimization Benchmark ===\n\n");
    
    // Run benchmark of field operations
    benchmark_field_operations();
    
    printf("\n=== Benchmark Complete ===\n");
    printf("For complete results, see the OPTIMIZATION_SUMMARY.md file\n");
    
    return 0;
}