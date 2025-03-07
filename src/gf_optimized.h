#ifndef GF_OPTIMIZED_H
#define GF_OPTIMIZED_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// Definition for the finite field GF(101)
#define MODULO_GF 101

// Field element type
typedef struct {
    uint8_t value;
} GF;

// Precomputed tables for all possible field operations
// These are initialized at program startup
static uint8_t GF_ADD_TABLE[MODULO_GF][MODULO_GF];
static uint8_t GF_SUB_TABLE[MODULO_GF][MODULO_GF];
static uint8_t GF_MUL_TABLE[MODULO_GF][MODULO_GF];
static uint8_t GF_INV_TABLE[MODULO_GF];
static uint8_t GF_NEG_TABLE[MODULO_GF];

// Memory pool for batch inversion
#define GF_POOL_SIZE 128
static GF* gf_pool[GF_POOL_SIZE] = {NULL};
static size_t gf_pool_sizes[GF_POOL_SIZE] = {0};
static bool gf_pool_used[GF_POOL_SIZE] = {false};
static size_t gf_pool_init = 0;

// Flag to check if tables are initialized
static bool GF_TABLES_INITIALIZED = false;

// Initialize all lookup tables
static void gf_init_tables() {
    if (GF_TABLES_INITIALIZED) return;
    
    // Addition table
    for (int i = 0; i < MODULO_GF; i++) {
        for (int j = 0; j < MODULO_GF; j++) {
            GF_ADD_TABLE[i][j] = (i + j) % MODULO_GF;
        }
    }
    
    // Subtraction table
    for (int i = 0; i < MODULO_GF; i++) {
        for (int j = 0; j < MODULO_GF; j++) {
            GF_SUB_TABLE[i][j] = (i - j + MODULO_GF) % MODULO_GF;
        }
    }
    
    // Multiplication table
    for (int i = 0; i < MODULO_GF; i++) {
        for (int j = 0; j < MODULO_GF; j++) {
            GF_MUL_TABLE[i][j] = (i * j) % MODULO_GF;
        }
    }
    
    // Negation table
    for (int i = 0; i < MODULO_GF; i++) {
        GF_NEG_TABLE[i] = (MODULO_GF - i) % MODULO_GF;
    }
    
    // Inverse table using Fermat's Little Theorem
    GF_INV_TABLE[0] = 0; // Special case: 0 has no inverse
    for (int i = 1; i < MODULO_GF; i++) {
        // Calculate i^(p-2) mod p using exponentiation by squaring
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
    
    GF_TABLES_INITIALIZED = true;
}

// Constructor that ensures tables are initialized
static inline GF gf_new(int64_t value) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    
    int64_t tmp = value % MODULO_GF;
    if (tmp < 0) tmp += MODULO_GF;
    
    GF result = { (uint8_t)tmp };
    return result;
}

// Alias for gf_new
static inline GF f101(int64_t n) {
    return gf_new(n);
}

// Constant field elements
static inline GF gf_zero() {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { 0 };
    return result;
}

static inline GF gf_one() {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { 1 };
    return result;
}

// Field validation
static inline bool gf_in_field(GF a) {
    return a.value < MODULO_GF;
}

// Field equality
static inline bool gf_equal(GF a, GF b) {
    return a.value == b.value;
}

// Field addition using lookup table
static inline GF gf_add(GF a, GF b) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { GF_ADD_TABLE[a.value][b.value] };
    return result;
}

// Field subtraction using lookup table
static inline GF gf_sub(GF a, GF b) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { GF_SUB_TABLE[a.value][b.value] };
    return result;
}

// Field multiplication using lookup table
static inline GF gf_mul(GF a, GF b) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { GF_MUL_TABLE[a.value][b.value] };
    return result;
}

// Field negation using lookup table
static inline GF gf_neg(GF a) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { GF_NEG_TABLE[a.value] };
    return result;
}

// Field inversion using lookup table
static inline GF gf_inv(GF a) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    GF result = { GF_INV_TABLE[a.value] };
    return result;
}

// Field division
static inline GF gf_div(GF a, GF b) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    if (b.value == 0) {
        fprintf(stderr, "Division by zero in finite field\n");
        exit(EXIT_FAILURE);
    }
    return gf_mul(a, gf_inv(b));
}

// Exponentiation with lookup-based squaring
static inline GF gf_pow(GF base, uint64_t exp) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    
    if (exp == 0) return gf_one();
    if (exp == 1) return base;
    
    GF result = gf_one();
    
    while (exp > 0) {
        if (exp & 1) {
            result = gf_mul(result, base);
        }
        base = gf_mul(base, base);
        exp >>= 1;
    }
    
    return result;
}

// Initialize the GF memory pool
static void gf_init_pool() {
    if (gf_pool_init) return;
    memset(gf_pool, 0, sizeof(gf_pool));
    memset(gf_pool_sizes, 0, sizeof(gf_pool_sizes));
    memset(gf_pool_used, 0, sizeof(gf_pool_used));
    gf_pool_init = 1;
}

// Get memory from the GF pool for batch operations
static GF* gf_get_memory(size_t size) {
    if (!gf_pool_init) gf_init_pool();
    
    // Try to find a suitable buffer in the pool
    for (int i = 0; i < GF_POOL_SIZE; i++) {
        if (!gf_pool_used[i] && gf_pool[i] && gf_pool_sizes[i] >= size) {
            gf_pool_used[i] = true;
            return gf_pool[i];
        }
    }
    
    // Allocate a new buffer with alignment (for SIMD)
    void* ptr = malloc(size * sizeof(GF));
    if (!ptr) {
        fprintf(stderr, "Memory allocation failed in gf_get_memory\n");
        exit(EXIT_FAILURE);
    }
    
    // Find an empty slot in the pool
    for (int i = 0; i < GF_POOL_SIZE; i++) {
        if (!gf_pool[i]) {
            gf_pool[i] = (GF*)ptr;
            gf_pool_sizes[i] = size;
            gf_pool_used[i] = true;
            return gf_pool[i];
        }
    }
    
    // If the pool is full, just return the new memory
    return (GF*)ptr;
}

// Return memory to the GF pool
static void gf_return_memory(GF* ptr) {
    if (!gf_pool_init || !ptr) return;
    
    for (int i = 0; i < GF_POOL_SIZE; i++) {
        if (gf_pool[i] == ptr) {
            gf_pool_used[i] = false;
            return;
        }
    }
    
    // If not found in the pool, add it if there's space
    for (int i = 0; i < GF_POOL_SIZE; i++) {
        if (!gf_pool[i]) {
            gf_pool[i] = ptr;
            gf_pool_used[i] = false;
            gf_pool_sizes[i] = 0;  // Size is unknown
            return;
        }
    }
    
    // If the pool is full, free the memory
    free(ptr);
}

// Batch inversion using Montgomery's trick
// Efficiently computes inversions of multiple field elements at once
// input_array and output_array can be the same array for in-place inversion
static void gf_batch_inv(const GF* input_array, GF* output_array, size_t count) {
    if (!GF_TABLES_INITIALIZED) gf_init_tables();
    
    if (count == 0) return;
    
    // Handle special case of single element
    if (count == 1) {
        output_array[0] = gf_inv(input_array[0]);
        return;
    }
    
    // Allocate temporary arrays from pool
    GF* temps = gf_get_memory(count);
    
    // Phase 1: Compute products
    // temps[i] = input[0] * input[1] * ... * input[i]
    temps[0] = input_array[0];
    for (size_t i = 1; i < count; i++) {
        temps[i] = gf_mul(temps[i-1], input_array[i]);
    }
    
    // Compute inverse of the product of all elements
    GF all_inverse = gf_inv(temps[count-1]);
    
    // Phase 2: Unwind the chain to get individual inverses
    // Last element's inverse can be computed directly
    output_array[count-1] = gf_mul(temps[count-2], all_inverse);
    
    // Compute the rest of inverses in reverse order
    for (size_t i = count-2; i > 0; i--) {
        all_inverse = gf_mul(all_inverse, input_array[i+1]); // Multiply by the next input element
        output_array[i] = gf_mul(temps[i-1], all_inverse);
    }
    
    // First element's inverse
    output_array[0] = gf_mul(all_inverse, input_array[1]);
    
    // Return memory to pool
    gf_return_memory(temps);
}

// Batch inversion with result array allocation
static GF* gf_batch_inv_new(const GF* input_array, size_t count) {
    if (count == 0) return NULL;
    
    GF* result = gf_get_memory(count);
    gf_batch_inv(input_array, result, count);
    return result;
}

#endif // GF_OPTIMIZED_H