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

#endif // GF_OPTIMIZED_H