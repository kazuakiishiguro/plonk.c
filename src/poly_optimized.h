#ifndef POLY_OPTIMIZED_H
#define POLY_OPTIMIZED_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include "hf.h"
#include "gf_optimized.h"

// Define alignment for vectorization
#define POLY_ALIGN 32

// Cache line size optimization
#define CACHE_LINE_SIZE 64  // Common size on modern CPUs

// Polynomial structure with cache alignment
typedef struct {
    HF* coeffs;         // Coefficient array
    size_t len;         // Length of polynomial
    size_t capacity;    // Allocated capacity (for reuse)
} POLY;

// Memory pool for polynomial operations to avoid repeated malloc/free
#define POLY_POOL_SIZE 32
static HF* poly_pool[POLY_POOL_SIZE] = {NULL};
static size_t poly_pool_sizes[POLY_POOL_SIZE] = {0};
static bool poly_pool_used[POLY_POOL_SIZE] = {false};
static size_t poly_pool_init = 0;

// Initialize the memory pool
static void poly_init_pool() {
    if (poly_pool_init) return;
    memset(poly_pool, 0, sizeof(poly_pool));
    memset(poly_pool_sizes, 0, sizeof(poly_pool_sizes));
    memset(poly_pool_used, 0, sizeof(poly_pool_used));
    poly_pool_init = 1;
}

// Get memory from the pool or allocate new memory
static HF* poly_get_memory(size_t size) {
    if (!poly_pool_init) poly_init_pool();
    
    // Try to find a suitable buffer in the pool
    for (int i = 0; i < POLY_POOL_SIZE; i++) {
        if (!poly_pool_used[i] && poly_pool[i] && poly_pool_sizes[i] >= size) {
            poly_pool_used[i] = true;
            return poly_pool[i];
        }
    }
    
    // Allocate a new buffer with alignment
    size_t aligned_size = ((size * sizeof(HF) + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
    void* ptr;
    int result = posix_memalign(&ptr, CACHE_LINE_SIZE, aligned_size);
    if (result != 0 || !ptr) {
        fprintf(stderr, "Memory allocation failed in poly_get_memory\n");
        exit(EXIT_FAILURE);
    }
    
    // Find an empty slot in the pool
    for (int i = 0; i < POLY_POOL_SIZE; i++) {
        if (!poly_pool[i]) {
            poly_pool[i] = (HF*)ptr;
            poly_pool_sizes[i] = size;
            poly_pool_used[i] = true;
            return poly_pool[i];
        }
    }
    
    // If the pool is full, just return the new memory (will be leaked if not properly freed)
    return (HF*)ptr;
}

// Return memory to the pool
static void poly_return_memory(HF* ptr) {
    if (!poly_pool_init || !ptr) return;
    
    for (int i = 0; i < POLY_POOL_SIZE; i++) {
        if (poly_pool[i] == ptr) {
            poly_pool_used[i] = false;
            return;
        }
    }
    
    // If not found in the pool, add it if there's space
    for (int i = 0; i < POLY_POOL_SIZE; i++) {
        if (!poly_pool[i]) {
            poly_pool[i] = ptr;
            poly_pool_used[i] = false;
            poly_pool_sizes[i] = 0;  // Size is unknown
            return;
        }
    }
    
    // If the pool is full, free the memory
    free(ptr);
}

// Create a new polynomial with trailing zeros removed
static POLY poly_new_internal(const HF* coeffs, size_t len) {
    // Find actual length by trimming trailing zeros
    size_t actual_len = len;
    while (actual_len > 1 && hf_equal(coeffs[actual_len - 1], hf_zero())) {
        actual_len--;
    }
    
    // Allocate memory from pool
    POLY result;
    result.len = actual_len;
    result.capacity = actual_len;
    result.coeffs = poly_get_memory(actual_len);
    
    // Use memcpy for efficient data transfer
    memcpy(result.coeffs, coeffs, actual_len * sizeof(HF));
    
    return result;
}

// Public constructor
static inline POLY poly_new(const HF* coeffs, size_t len) {
    return poly_new_internal(coeffs, len);
}

// Zero polynomial
static inline POLY poly_zero() {
    HF zero = hf_zero();
    return poly_new(&zero, 1);
}

// One polynomial
static inline POLY poly_one() {
    HF one = hf_one();
    return poly_new(&one, 1);
}

// Check if polynomial is zero
static inline bool poly_is_zero(const POLY* polynomial) {
    if (polynomial->len == 0) return true;
    if (polynomial->len == 1 && hf_equal(polynomial->coeffs[0], hf_zero())) return true;
    return false;
}

// Add a scalar to a polynomial (modifies in place)
static inline void poly_add_hf_inplace(POLY* a, const HF b) {
    a->coeffs[0] = hf_add(a->coeffs[0], b);
}

// Add scalar and return new polynomial
static inline POLY poly_add_hf(const POLY* a, const HF b) {
    POLY result = poly_new(a->coeffs, a->len);
    poly_add_hf_inplace(&result, b);
    return result;
}

// Add two polynomials with vectorization-friendly loop structure
static inline POLY poly_add(const POLY* a, const POLY* b) {
    size_t max_len = (a->len > b->len) ? a->len : b->len;
    HF* coeffs = poly_get_memory(max_len);
    
    // Process the common part of both polynomials
    size_t common_len = (a->len < b->len) ? a->len : b->len;
    
    // Main vectorizable loop
    for (size_t i = 0; i < common_len; i++) {
        coeffs[i] = hf_add(a->coeffs[i], b->coeffs[i]);
    }
    
    // Copy the rest from the longer polynomial
    if (a->len > common_len) {
        memcpy(coeffs + common_len, a->coeffs + common_len, (a->len - common_len) * sizeof(HF));
    } else if (b->len > common_len) {
        memcpy(coeffs + common_len, b->coeffs + common_len, (b->len - common_len) * sizeof(HF));
    }
    
    POLY result;
    result.coeffs = coeffs;
    result.len = max_len;
    result.capacity = max_len;
    
    // Remove trailing zeros if needed
    while (result.len > 1 && hf_equal(result.coeffs[result.len - 1], hf_zero())) {
        result.len--;
    }
    
    return result;
}

// Subtract polynomials
static inline POLY poly_sub(const POLY* a, const POLY* b) {
    size_t max_len = (a->len > b->len) ? a->len : b->len;
    HF* coeffs = poly_get_memory(max_len);
    
    // Process common part
    size_t common_len = (a->len < b->len) ? a->len : b->len;
    for (size_t i = 0; i < common_len; i++) {
        coeffs[i] = hf_sub(a->coeffs[i], b->coeffs[i]);
    }
    
    // Copy rest from longer polynomial (with negation if from b)
    if (a->len > common_len) {
        memcpy(coeffs + common_len, a->coeffs + common_len, (a->len - common_len) * sizeof(HF));
    } else if (b->len > common_len) {
        // Need to negate coefficients from b
        for (size_t i = common_len; i < b->len; i++) {
            coeffs[i] = hf_neg(b->coeffs[i]);
        }
    }
    
    POLY result;
    result.coeffs = coeffs;
    result.len = max_len;
    result.capacity = max_len;
    
    // Remove trailing zeros
    while (result.len > 1 && hf_equal(result.coeffs[result.len - 1], hf_zero())) {
        result.len--;
    }
    
    return result;
}

// Multiply polynomials with cache-optimized blocked algorithm
static inline POLY poly_mul(const POLY* a, const POLY* b) {
    // Handle special cases
    if (poly_is_zero(a) || poly_is_zero(b)) return poly_zero();
    if (a->len == 1 && hf_equal(a->coeffs[0], hf_one())) return poly_new(b->coeffs, b->len);
    if (b->len == 1 && hf_equal(b->coeffs[0], hf_one())) return poly_new(a->coeffs, a->len);
    
    // Resulting polynomial has degree len_a + len_b - 1
    size_t result_len = a->len + b->len - 1;
    HF* coeffs = poly_get_memory(result_len);
    memset(coeffs, 0, result_len * sizeof(HF));  // Initialize to zero
    
    // Cache blocking parameters
    const size_t BLOCK_SIZE = 16;  // Tune based on your cache size
    
    // Blocked multiplication for better cache locality
    for (size_t i_block = 0; i_block < a->len; i_block += BLOCK_SIZE) {
        size_t i_end = (i_block + BLOCK_SIZE < a->len) ? i_block + BLOCK_SIZE : a->len;
        
        for (size_t j_block = 0; j_block < b->len; j_block += BLOCK_SIZE) {
            size_t j_end = (j_block + BLOCK_SIZE < b->len) ? j_block + BLOCK_SIZE : b->len;
            
            // Multiply blocks
            for (size_t i = i_block; i < i_end; i++) {
                const HF a_coeff = a->coeffs[i];
                if (hf_equal(a_coeff, hf_zero())) continue;
                
                for (size_t j = j_block; j < j_end; j++) {
                    const HF b_coeff = b->coeffs[j];
                    if (hf_equal(b_coeff, hf_zero())) continue;
                    
                    const HF product = hf_mul(a_coeff, b_coeff);
                    coeffs[i + j] = hf_add(coeffs[i + j], product);
                }
            }
        }
    }
    
    POLY result;
    result.coeffs = coeffs;
    result.len = result_len;
    result.capacity = result_len;
    
    // Remove trailing zeros if needed
    while (result.len > 1 && hf_equal(result.coeffs[result.len - 1], hf_zero())) {
        result.len--;
    }
    
    return result;
}

// Scale a polynomial by a constant
static inline POLY poly_scale(const POLY* p, HF scalar) {
    if (hf_equal(scalar, hf_zero())) return poly_zero();
    if (hf_equal(scalar, hf_one())) return poly_new(p->coeffs, p->len);
    
    HF* coeffs = poly_get_memory(p->len);
    
    // Vectorizable loop
    for (size_t i = 0; i < p->len; i++) {
        coeffs[i] = hf_mul(p->coeffs[i], scalar);
    }
    
    POLY result;
    result.coeffs = coeffs;
    result.len = p->len;
    result.capacity = p->len;
    
    // Remove trailing zeros if needed
    while (result.len > 1 && hf_equal(result.coeffs[result.len - 1], hf_zero())) {
        result.len--;
    }
    
    return result;
}

// Evaluate polynomial using Horner's method with vectorization-friendly structure
static inline HF poly_eval(const POLY* p, HF x) {
    if (p->len == 0) return hf_zero();
    if (p->len == 1) return p->coeffs[0];
    if (hf_equal(x, hf_zero())) return p->coeffs[0];
    if (hf_equal(x, hf_one())) {
        // Sum all coefficients
        HF sum = hf_zero();
        for (size_t i = 0; i < p->len; i++) {
            sum = hf_add(sum, p->coeffs[i]);
        }
        return sum;
    }
    
    // Standard Horner's method with unrolling
    HF result = p->coeffs[p->len - 1];
    
    // Process 4 coefficients at a time for better instruction parallelism
    size_t i = p->len - 2;
    for (; i >= 3; i -= 4) {
        result = hf_mul(result, x);
        result = hf_add(result, p->coeffs[i]);
        
        result = hf_mul(result, x);
        result = hf_add(result, p->coeffs[i-1]);
        
        result = hf_mul(result, x);
        result = hf_add(result, p->coeffs[i-2]);
        
        result = hf_mul(result, x);
        result = hf_add(result, p->coeffs[i-3]);
    }
    
    // Handle remaining coefficients
    for (; i < p->len; i--) {
        result = hf_mul(result, x);
        result = hf_add(result, p->coeffs[i]);
    }
    
    return result;
}

// Free a polynomial and return memory to pool
static inline void poly_free(POLY* p) {
    if (p && p->coeffs) {
        poly_return_memory(p->coeffs);
        p->coeffs = NULL;
        p->len = 0;
        p->capacity = 0;
    }
}

#endif // POLY_OPTIMIZED_H