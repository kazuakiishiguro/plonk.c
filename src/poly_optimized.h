#ifndef POLY_OPTIMIZED_H
#define POLY_OPTIMIZED_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "hf.h"
#include "gf_optimized.h"

// Define alignment for vectorization
#define POLY_ALIGN 32

// Cache line size optimization
#define CACHE_LINE_SIZE 64  // Common size on modern CPUs

// FFT threshold - use FFT for polynomials larger than this
#define FFT_THRESHOLD 64

// Use power of 2 sizes for FFT efficiency
#define FFT_MIN_SIZE 64   // Minimum FFT size (must be power of 2)

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

// Find the next power of 2 greater than or equal to n
static inline size_t next_power_of_two(size_t n) {
    size_t power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

// Fast Fourier Transform for polynomial multiplication
// This is a recursive implementation of the Cooley-Tukey algorithm
static void poly_fft(HF* in, HF* out, size_t n, HF omega, bool inverse) {
    if (n == 1) {
        out[0] = in[0];
        return;
    }
    
    // Split into even and odd indices
    HF* even = poly_get_memory(n/2);
    HF* odd = poly_get_memory(n/2);
    
    for (size_t i = 0; i < n/2; i++) {
        even[i] = in[2*i];
        odd[i] = in[2*i + 1];
    }
    
    // Allocate space for results of recursive calls
    HF* even_fft = poly_get_memory(n/2);
    HF* odd_fft = poly_get_memory(n/2);
    
    // Recursively compute FFT of even and odd parts
    HF omega_squared = hf_mul(omega, omega);
    poly_fft(even, even_fft, n/2, omega_squared, inverse);
    poly_fft(odd, odd_fft, n/2, omega_squared, inverse);
    
    // Combine results
    HF omega_pow = hf_one();
    for (size_t i = 0; i < n/2; i++) {
        HF t = hf_mul(omega_pow, odd_fft[i]);
        out[i] = hf_add(even_fft[i], t);
        out[i + n/2] = hf_sub(even_fft[i], t);
        omega_pow = hf_mul(omega_pow, omega);
    }
    
    // Clean up temporary memory
    poly_return_memory(even);
    poly_return_memory(odd);
    poly_return_memory(even_fft);
    poly_return_memory(odd_fft);
}

// Compute primitive root of unity for FFT
static HF get_omega(size_t n, bool inverse) {
    // For GF(101), we can use 2 as a primitive root
    // We need an element of order n
    size_t order = MODULO_GF - 1; // order of the multiplicative group
    size_t k = order / n;         // we need an element of order n, which is (p-1)/n
    
    HF omega = hf_pow(hf_new(2), k); // 2^k is a primitive n-th root of unity
    
    if (inverse) {
        // For inverse FFT, we need the inverse of omega
        omega = hf_inv(omega);
    }
    
    return omega;
}

// Multiply polynomials using FFT for large polynomials, or standard algorithm for small ones
static inline POLY poly_mul(const POLY* a, const POLY* b) {
    // Handle special cases
    if (poly_is_zero(a) || poly_is_zero(b)) return poly_zero();
    if (a->len == 1 && hf_equal(a->coeffs[0], hf_one())) return poly_new(b->coeffs, b->len);
    if (b->len == 1 && hf_equal(b->coeffs[0], hf_one())) return poly_new(a->coeffs, a->len);
    
    // Resulting polynomial has degree len_a + len_b - 1
    size_t result_len = a->len + b->len - 1;
    
    // Use classic algorithm for small polynomials
    if (result_len < FFT_THRESHOLD) {
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
    
    // For larger polynomials, use FFT-based multiplication
    // Find power of 2 size for FFT that can hold both polynomials
    size_t fft_size = next_power_of_two(result_len);
    
    // Prepare input arrays with zero padding
    HF* a_padded = poly_get_memory(fft_size);
    HF* b_padded = poly_get_memory(fft_size);
    
    // Zero-fill
    memset(a_padded, 0, fft_size * sizeof(HF));
    memset(b_padded, 0, fft_size * sizeof(HF));
    
    // Copy coefficients
    memcpy(a_padded, a->coeffs, a->len * sizeof(HF));
    memcpy(b_padded, b->coeffs, b->len * sizeof(HF));
    
    // Allocate FFT result arrays
    HF* a_fft = poly_get_memory(fft_size);
    HF* b_fft = poly_get_memory(fft_size);
    
    // Compute FFT of input polynomials
    HF omega = get_omega(fft_size, false);
    poly_fft(a_padded, a_fft, fft_size, omega, false);
    poly_fft(b_padded, b_fft, fft_size, omega, false);
    
    // Point-wise multiplication of FFT results
    HF* c_fft = poly_get_memory(fft_size);
    for (size_t i = 0; i < fft_size; i++) {
        c_fft[i] = hf_mul(a_fft[i], b_fft[i]);
    }
    
    // Inverse FFT to get result coefficients
    HF* c_coeffs = poly_get_memory(fft_size);
    poly_fft(c_fft, c_coeffs, fft_size, get_omega(fft_size, true), true);
    
    // Normalize by dividing by fft_size (we need to use batch inversion here)
    HF n_inv = hf_inv(hf_new(fft_size % MODULO_GF));
    for (size_t i = 0; i < fft_size; i++) {
        c_coeffs[i] = hf_mul(c_coeffs[i], n_inv);
    }
    
    // Create result polynomial with the correct size
    POLY result;
    result.len = result_len;
    result.capacity = fft_size;
    result.coeffs = c_coeffs;
    
    // Remove trailing zeros if needed
    while (result.len > 1 && hf_equal(result.coeffs[result.len - 1], hf_zero())) {
        result.len--;
    }
    
    // Free temporary arrays
    poly_return_memory(a_padded);
    poly_return_memory(b_padded);
    poly_return_memory(a_fft);
    poly_return_memory(b_fft);
    poly_return_memory(c_fft);
    
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

// Optimized vanishing polynomial using efficient memory management and batch operations
static POLY poly_z(const HF *points, size_t len) {
    if (len == 0) return poly_one();
    if (len == 1) {
        HF coeffs[] = { hf_neg(points[0]), hf_one() };
        return poly_new(coeffs, 2);
    }
    
    // For small lengths, use the standard approach
    if (len < 8) {
        // Construct Π (x - points[i])
        POLY acc = poly_one();
        for (size_t i = 0; i < len; i++) {
            HF coeffs[] = { hf_neg(points[i]), hf_one() };
            POLY term = poly_new(coeffs, 2);
            POLY temp = poly_mul(&acc, &term);
            poly_free(&acc);
            poly_free(&term);
            acc = temp;
        }
        return acc;
    }
    
    // For larger lengths, use a divide-and-conquer approach
    // This reduces the number of polynomial multiplications
    size_t half = len / 2;
    
    // Compute left and right halves recursively
    POLY left = poly_z(points, half);
    POLY right = poly_z(points + half, len - half);
    
    // Multiply the results
    POLY result = poly_mul(&left, &right);
    
    // Clean up
    poly_free(&left);
    poly_free(&right);
    
    return result;
}

// Optimized Lagrange interpolation with batch inversion
static POLY poly_lagrange(const HF *x_points, const HF *y_points, size_t len) {
    // Special cases
    if (len == 0) return poly_zero();
    if (len == 1) {
        POLY result = poly_zero();
        result.coeffs[0] = y_points[0];
        return result;
    }
    
    // Compute all differences between x points
    // (x_j - x_i) for all i,j where i != j
    size_t num_diffs = len * (len - 1);
    HF* differences = poly_get_memory(num_diffs);
    size_t* diff_indices = (size_t*)malloc(num_diffs * 2 * sizeof(size_t));
    
    if (!diff_indices) {
        fprintf(stderr, "Memory allocation failed in poly_lagrange\n");
        exit(EXIT_FAILURE);
    }
    
    size_t diff_count = 0;
    for (size_t j = 0; j < len; j++) {
        for (size_t i = 0; i < len; i++) {
            if (i != j) {
                // Store the difference (x_j - x_i)
                differences[diff_count] = hf_sub(x_points[j], x_points[i]);
                
                // Store indices for lookup
                diff_indices[diff_count*2] = j;
                diff_indices[diff_count*2 + 1] = i;
                diff_count++;
            }
        }
    }
    
    // Compute all inverses at once using batch inversion
    HF* inverses = poly_get_memory(num_diffs);
    // Convert from HF to GF for batch inversion
    GF* gf_diffs = (GF*)malloc(num_diffs * sizeof(GF));
    GF* gf_invs = (GF*)malloc(num_diffs * sizeof(GF));
    
    if (!gf_diffs || !gf_invs) {
        fprintf(stderr, "Memory allocation failed in poly_lagrange\n");
        exit(EXIT_FAILURE);
    }
    
    // Prepare inputs for batch inversion
    for (size_t i = 0; i < num_diffs; i++) {
        gf_diffs[i] = differences[i].value;
    }
    
    // Batch invert all differences
    gf_batch_inv(gf_diffs, gf_invs, num_diffs);
    
    // Convert back to HF format
    for (size_t i = 0; i < num_diffs; i++) {
        inverses[i].value = gf_invs[i].value;
    }
    
    free(gf_diffs);
    free(gf_invs);
    
    // Now build the Lagrange polynomial with all inverses available
    POLY l = poly_zero();
    
    for (size_t j = 0; j < len; j++) {
        POLY l_j = poly_one();
        
        for (size_t i = 0; i < len; i++) {
            if (i != j) {
                // Find the precomputed inverse
                size_t inv_idx = 0;
                for (size_t k = 0; k < num_diffs; k++) {
                    if (diff_indices[k*2] == j && diff_indices[k*2 + 1] == i) {
                        inv_idx = k;
                        break;
                    }
                }
                
                // Get cached inverse
                HF denom_inv = inverses[inv_idx];
                
                HF neg_cxi = hf_neg(hf_mul(denom_inv, x_points[i]));
                HF coeffs[] = {neg_cxi, denom_inv};
                POLY term = poly_new(coeffs, 2);
                POLY temp = poly_mul(&l_j, &term);
                poly_free(&l_j);
                poly_free(&term);
                l_j = temp;
            }
        }
        
        POLY scaled_l_j = poly_scale(&l_j, y_points[j]);
        POLY temp_l = poly_add(&l, &scaled_l_j);
        poly_free(&l);
        poly_free(&l_j);
        poly_free(&scaled_l_j);
        l = temp_l;
    }
    
    // Clean up
    poly_return_memory(differences);
    poly_return_memory(inverses);
    free(diff_indices);
    
    return l;
}

#endif // POLY_OPTIMIZED_H