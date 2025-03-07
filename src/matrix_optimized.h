#ifndef MATRIX_OPTIMIZED_H
#define MATRIX_OPTIMIZED_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hf.h"
#include "gf_optimized.h"

// Define optimization parameters
#define MATRIX_ALIGN 64
#define MATRIX_BLOCK_SIZE 16

// Matrix structure with cache-friendly array layout
typedef struct {
    size_t m;           // Number of rows
    size_t n;           // Number of columns
    size_t stride;      // Stride for aligned access (may be > n)
    HF* v;              // Flat array of field elements
} MATRIX;

// Memory pool for matrices to reduce allocations
#define MATRIX_POOL_SIZE 16
static HF* matrix_pool[MATRIX_POOL_SIZE] = {NULL};
static size_t matrix_pool_sizes[MATRIX_POOL_SIZE] = {0};
static bool matrix_pool_used[MATRIX_POOL_SIZE] = {false};
static size_t matrix_pool_init = 0;

// Initialize the matrix memory pool
static void matrix_init_pool() {
    if (matrix_pool_init) return;
    memset(matrix_pool, 0, sizeof(matrix_pool));
    memset(matrix_pool_sizes, 0, sizeof(matrix_pool_sizes));
    memset(matrix_pool_used, 0, sizeof(matrix_pool_used));
    matrix_pool_init = 1;
}

// Get aligned memory from the pool or allocate new memory
static HF* matrix_get_memory(size_t size) {
    if (!matrix_pool_init) matrix_init_pool();
    
    // Try to find suitable memory in the pool
    for (int i = 0; i < MATRIX_POOL_SIZE; i++) {
        if (!matrix_pool_used[i] && matrix_pool[i] && matrix_pool_sizes[i] >= size) {
            matrix_pool_used[i] = true;
            return matrix_pool[i];
        }
    }
    
    // Allocate new aligned memory
    void* ptr;
    int result = posix_memalign(&ptr, MATRIX_ALIGN, size * sizeof(HF));
    if (result != 0 || !ptr) {
        fprintf(stderr, "Memory allocation failed in matrix_get_memory\n");
        exit(EXIT_FAILURE);
    }
    
    // Store in pool if possible
    for (int i = 0; i < MATRIX_POOL_SIZE; i++) {
        if (!matrix_pool[i]) {
            matrix_pool[i] = (HF*)ptr;
            matrix_pool_sizes[i] = size;
            matrix_pool_used[i] = true;
            return matrix_pool[i];
        }
    }
    
    // If pool is full, just return the memory
    return (HF*)ptr;
}

// Return memory to the pool
static void matrix_return_memory(HF* ptr) {
    if (!matrix_pool_init || !ptr) return;
    
    for (int i = 0; i < MATRIX_POOL_SIZE; i++) {
        if (matrix_pool[i] == ptr) {
            matrix_pool_used[i] = false;
            return;
        }
    }
    
    // If not in pool, try to add it
    for (int i = 0; i < MATRIX_POOL_SIZE; i++) {
        if (!matrix_pool[i]) {
            matrix_pool[i] = ptr;
            matrix_pool_used[i] = false;
            return;
        }
    }
    
    // If pool is full, free the memory
    free(ptr);
}

// Create a zero matrix with optimized memory layout
static MATRIX matrix_zero(size_t m, size_t n) {
    // Calculate stride for alignment
    size_t stride = ((n + MATRIX_ALIGN - 1) / MATRIX_ALIGN) * MATRIX_ALIGN;
    
    // Allocate aligned memory
    size_t size = m * stride;
    HF* v = matrix_get_memory(size);
    
    // Initialize to zero
    memset(v, 0, size * sizeof(HF));
    
    MATRIX result;
    result.m = m;
    result.n = n;
    result.stride = stride;
    result.v = v;
    
    return result;
}

// Create a new matrix from an array of values
static MATRIX matrix_new(HF* values, size_t m, size_t n) {
    // Create a zero matrix with proper alignment
    MATRIX result = matrix_zero(m, n);
    
    // Copy values with potential stride adjustment
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            result.v[i * result.stride + j] = values[i * n + j];
        }
    }
    
    return result;
}

// Get matrix element with bounds checking
static inline HF matrix_get(const MATRIX* matrix, size_t row, size_t col) {
    if (row >= matrix->m || col >= matrix->n) {
        fprintf(stderr, "Index out of bounds in matrix_get: [%zu,%zu] in [%zu,%zu]\n", 
                row, col, matrix->m, matrix->n);
        exit(EXIT_FAILURE);
    }
    return matrix->v[row * matrix->stride + col];
}

// Set matrix element with bounds checking
static inline void matrix_set(MATRIX* matrix, size_t row, size_t col, HF value) {
    if (row >= matrix->m || col >= matrix->n) {
        fprintf(stderr, "Index out of bounds in matrix_set\n");
        exit(EXIT_FAILURE);
    }
    matrix->v[row * matrix->stride + col] = value;
}

// Free a matrix and return memory to pool
static inline void matrix_free(MATRIX* matrix) {
    if (matrix->v) {
        matrix_return_memory(matrix->v);
        matrix->v = NULL;
    }
    matrix->m = 0;
    matrix->n = 0;
    matrix->stride = 0;
}

// Add two matrices with vectorizable loops
static MATRIX matrix_add(const MATRIX* a, const MATRIX* b) {
    if (a->m != b->m || a->n != b->n) {
        fprintf(stderr, "Matrix dimensions must match for addition\n");
        exit(EXIT_FAILURE);
    }
    
    // Create result matrix
    MATRIX result = matrix_zero(a->m, a->n);
    
    // Add matrices with stride adjustment for better cache behavior
    for (size_t i = 0; i < a->m; i++) {
        for (size_t j = 0; j < a->n; j++) {
            size_t a_idx = i * a->stride + j;
            size_t b_idx = i * b->stride + j;
            size_t res_idx = i * result.stride + j;
            result.v[res_idx] = hf_add(a->v[a_idx], b->v[b_idx]);
        }
    }
    
    return result;
}

// Multiply matrices with cache-blocking optimization
static MATRIX matrix_mul(const MATRIX* a, const MATRIX* b) {
    if (a->n != b->m) {
        fprintf(stderr, "Matrix dimensions incompatible for multiplication\n");
        exit(EXIT_FAILURE);
    }
    
    // Create result matrix
    MATRIX result = matrix_zero(a->m, b->n);
    
    // Cache blocking parameters
    const size_t block_size = MATRIX_BLOCK_SIZE;
    
    // Blocked matrix multiplication
    for (size_t i_block = 0; i_block < a->m; i_block += block_size) {
        size_t i_end = (i_block + block_size < a->m) ? i_block + block_size : a->m;
        
        for (size_t j_block = 0; j_block < b->n; j_block += block_size) {
            size_t j_end = (j_block + block_size < b->n) ? j_block + block_size : b->n;
            
            for (size_t k_block = 0; k_block < a->n; k_block += block_size) {
                size_t k_end = (k_block + block_size < a->n) ? k_block + block_size : a->n;
                
                // Process block
                for (size_t i = i_block; i < i_end; i++) {
                    for (size_t j = j_block; j < j_end; j++) {
                        HF sum = hf_zero();
                        
                        // Inner loop with direct memory access
                        for (size_t k = k_block; k < k_end; k++) {
                            HF a_val = a->v[i * a->stride + k];
                            HF b_val = b->v[k * b->stride + j];
                            sum = hf_add(sum, hf_mul(a_val, b_val));
                        }
                        
                        // Accumulate into result
                        size_t res_idx = i * result.stride + j;
                        result.v[res_idx] = hf_add(result.v[res_idx], sum);
                    }
                }
            }
        }
    }
    
    return result;
}

// Efficient matrix inversion with Gauss-Jordan elimination
static MATRIX matrix_inv(const MATRIX* matrix) {
    if (matrix->m != matrix->n) {
        fprintf(stderr, "Only square matrices can be inverted\n");
        exit(EXIT_FAILURE);
    }
    
    size_t n = matrix->n;
    
    // Create augmented matrix [A|I]
    MATRIX aug = matrix_zero(n, 2 * n);
    
    // Copy input matrix to left side
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            matrix_set(&aug, i, j, matrix_get(matrix, i, j));
        }
        // Set identity matrix on right side
        matrix_set(&aug, i, i + n, hf_one());
    }
    
    // Perform Gauss-Jordan elimination
    for (size_t i = 0; i < n; i++) {
        // Find pivot
        size_t pivot_row = i;
        for (size_t j = i + 1; j < n; j++) {
            if (!hf_equal(matrix_get(&aug, j, i), hf_zero())) {
                pivot_row = j;
                break;
            }
        }
        
        // Swap rows if needed
        if (pivot_row != i) {
            for (size_t j = 0; j < 2 * n; j++) {
                HF temp = matrix_get(&aug, i, j);
                matrix_set(&aug, i, j, matrix_get(&aug, pivot_row, j));
                matrix_set(&aug, pivot_row, j, temp);
            }
        }
        
        // Scale pivot row
        HF pivot = matrix_get(&aug, i, i);
        if (hf_equal(pivot, hf_zero())) {
            fprintf(stderr, "Matrix is singular and cannot be inverted\n");
            exit(EXIT_FAILURE);
        }
        
        HF pivot_inv = hf_inv(pivot);
        for (size_t j = 0; j < 2 * n; j++) {
            matrix_set(&aug, i, j, hf_mul(matrix_get(&aug, i, j), pivot_inv));
        }
        
        // Eliminate other rows
        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;
            
            HF factor = matrix_get(&aug, j, i);
            for (size_t k = 0; k < 2 * n; k++) {
                HF subtrahend = hf_mul(matrix_get(&aug, i, k), factor);
                matrix_set(&aug, j, k, hf_sub(matrix_get(&aug, j, k), subtrahend));
            }
        }
    }
    
    // Extract inverse from right side of augmented matrix
    MATRIX inv = matrix_zero(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            matrix_set(&inv, i, j, matrix_get(&aug, i, j + n));
        }
    }
    
    // Free augmented matrix
    matrix_free(&aug);
    
    return inv;
}

#endif // MATRIX_OPTIMIZED_H