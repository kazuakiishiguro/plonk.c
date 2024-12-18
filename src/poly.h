#ifndef POLY_H
#define POLY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include "hf.h"

typedef struct {
  HF *coeffs;
  size_t len;
} POLY;

/**
 * Create a polynomial and automatically remove trailing zeros.
 * coeffs: Input array of coefficients.
 * len: Length of the input array.
 */
static POLY poly_new_internal(const HF *coeffs, size_t len) {
  // Trim trailing zeros
  while (len > 1 && hf_equal(coeffs[len - 1], hf_zero())) {
    len--;
  }

  POLY result;
  result.len = len;
  result.coeffs = (HF *)malloc(len * sizeof(HF));
  if (!result.coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_new_internal\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < len; i++) {
    result.coeffs[i] = coeffs[i];
  }
  return result;
}

// Public constructor-like function
static inline POLY poly_new(const HF *coeffs, size_t len) {
  return poly_new_internal(coeffs, len);
}

static inline POLY poly_zero() {
  HF zero = hf_zero();
  return poly_new(&zero, 1);
}

static inline POLY poly_one() {
  HF one = hf_one();
  return poly_new(&one, 1);
}

static inline bool poly_is_zero(const POLY *polynomial) {
  if (polynomial->len == 0) {
    return true;
  }
  for (size_t i = 0; i < polynomial->len; i++) {
    if (!hf_equal(polynomial->coeffs[i], hf_zero()))
      return false;
  }
  return true;
}

// Add a scalar HF value to the polynomial (to the constant term)
POLY poly_add_hf(POLY *a, const HF b) {
  a->coeffs[0] = hf_add(a->coeffs[0], b);
  return *a;
}

POLY poly_add(const POLY *a, const POLY *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  HF *coeffs = (HF *)malloc(max_len * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_add\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    HF a_coeff = (i < a->len) ? a->coeffs[i] : hf_new(0);
    HF b_coeff = (i < b->len) ? b->coeffs[i] : hf_new(0);
    coeffs[i] = hf_add(a_coeff, b_coeff);
  }
  POLY result = poly_new(coeffs, max_len);
  free(coeffs);
  return result;
}

POLY poly_sub(const POLY *a, const POLY *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  HF *coeffs = (HF *)malloc(max_len * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_sub\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    HF a_coeff = (i < a->len) ? a->coeffs[i] : hf_new(0);
    HF b_coeff = (i < b->len) ? b->coeffs[i] : hf_new(0);
    coeffs[i] = hf_sub(a_coeff, b_coeff);
  }
  POLY p = poly_new(coeffs, max_len);
  free(coeffs);
  return p;
}

POLY poly_mul(const POLY *a, const POLY *b) {
  size_t result_len = a->len + b->len - 1;
  HF *coeffs = (HF *)calloc(result_len, sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_mul\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < a->len; i++) {
    for (size_t j = 0; j < b->len; j++) {
      HF product = hf_mul(a->coeffs[i], b->coeffs[j]);
      coeffs[i + j] = hf_add(coeffs[i + j], product);
    }
  }
  POLY reuslt = poly_new(coeffs, result_len);
  free(coeffs);
  return reuslt;
}

void poly_divide(const POLY *num, const POLY *den, POLY *quot, POLY *rem) {
  if (poly_is_zero(den)) {
    fprintf(stderr, "Division by zero polynomial in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  size_t num_len = num->len;
  size_t den_len = den->len;

  HF *quot_coeffs = (HF *)calloc(num_len, sizeof(HF));
  HF *rem_coeffs  = (HF *)calloc(num_len, sizeof(HF));
  if (!quot_coeffs || !rem_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // Copy numerator -> remainder
  for (size_t i = 0; i < num_len; i++) {
    rem_coeffs[i] = num->coeffs[i];
  }

  HF den_lead_inv = hf_inv(den->coeffs[den_len - 1]);

  for (ssize_t i = (ssize_t)num_len - 1; i >= (ssize_t)(den_len - 1); i--) {
    HF coeff = hf_mul(rem_coeffs[i], den_lead_inv);
    quot_coeffs[i - (den_len - 1)] = coeff;

    for (ssize_t j = 0; j < (ssize_t)den_len; j++) {
      rem_coeffs[i - j] = hf_sub(rem_coeffs[i - j], hf_mul(coeff, den->coeffs[den_len - 1 - j]));
    }
  }

  // Trim quotient
  size_t quot_len = num_len >= den_len ? num_len - den_len + 1 : 1;
  while (quot_len > 1 && hf_equal(quot_coeffs[quot_len - 1], hf_zero())) {
    quot_len--;
  }

  // Trim remainder
  size_t rem_len = den_len - 1;
  if (rem_len > num_len) {
    // If denominator is longer than numerator, remainder is actually just numerator.
    rem_len = num_len;
  }
  while (rem_len > 1 && hf_equal(rem_coeffs[rem_len - 1], hf_zero())) {
    rem_len--;
  }

  *quot = poly_new(quot_coeffs, quot_len);
  *rem  = poly_new(rem_coeffs,  rem_len);

  free(quot_coeffs);
  free(rem_coeffs);
}

POLY poly_scale(const POLY *p, HF scalar) {
  if (hf_equal(scalar, hf_zero())) {
    return poly_zero();
  }

  HF *coeffs = (HF *)malloc(p->len * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_scale\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < p->len; i++) {
    coeffs[i] = hf_mul(p->coeffs[i], scalar);
  }

  POLY result = poly_new(coeffs, p->len);
  free(coeffs);
  return result;
}

POLY poly_shift(const POLY *p, size_t shift) {
  if (poly_is_zero(p)) return poly_zero();

  size_t new_len = p->len + shift;
  HF *new_coeffs = (HF *)calloc(new_len, sizeof(HF));
  if (!new_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_shift\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < p->len; i++) {
    new_coeffs[i + shift] = p->coeffs[i];
  }

  POLY result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

POLY poly_slice(const POLY *p, size_t start, size_t end) {
  if (start >= end || end > p->len) {
    fprintf(stderr, "Invalid slice indices in poly_slice\n");
    exit(EXIT_FAILURE);
  }

  size_t new_len = end - start;
  HF *new_coeffs = (HF *)malloc(new_len * sizeof(HF));
  if (!new_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_slice\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < new_len; i++) {
    new_coeffs[i] = p->coeffs[start + i];
  }

  POLY result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

POLY poly_negate(const POLY *p) {
  HF *neg_coeffs = (HF *)malloc(p->len * sizeof(HF));
  if (!neg_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_negate\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < p->len; i++) {
    neg_coeffs[i] = hf_neg(p->coeffs[i]);
  }

  POLY result = poly_new(neg_coeffs, p->len);
  free(neg_coeffs);
  return result;
}

static inline void poly_free(POLY *p) {
  if (p->coeffs) {
    free(p->coeffs);
    p->coeffs = NULL;
  }
  p->len = 0;
}

// Evaluate polynomial using Horner's method: O(n) with minimal overhead
HF poly_eval(const POLY *p, HF x) {
  HF y = hf_zero();
  for (ssize_t i = (ssize_t)p->len - 1; i >= 0; i--) {
    y = hf_mul(y, x);
    y = hf_add(y, p->coeffs[i]);
  }
  return y;
}

POLY poly_z(const HF *points, size_t len) {
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

POLY poly_lagrange(const HF *x_points, const HF *y_points, size_t len) {
  // Lagrange interpolation:
  // L(x) = Σ [y_j * Π((x - x_i)/(x_j - x_i)) for i != j]
  POLY l = poly_zero();
  for (size_t j = 0; j < len; j++) {
    POLY l_j = poly_one();
    for (size_t i = 0; i < len; i++) {
      if (i != j) {
        HF denom = hf_sub(x_points[j], x_points[i]);
        HF denom_inv = hf_inv(denom);
        if (hf_equal(denom_inv, hf_zero())) {
          fprintf(stderr, "Error: Lagrange polynomial x points must be unique\n");
          exit(EXIT_FAILURE);
        }

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
  return l;
}

#endif // POLY_H
