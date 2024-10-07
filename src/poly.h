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

void poly_normalize(POLY *polynomial) {
  while (polynomial->len > 1 && hf_equal(polynomial->coeffs[polynomial->len - 1], hf_new(0))) {
    polynomial->len--;
  }
  polynomial->coeffs = realloc(polynomial->coeffs, polynomial->len * sizeof(HF));
};

POLY poly_new(HF *coeffs, size_t len) {
  POLY result;
  result.len = len;
  result.coeffs = (HF *)malloc(len * sizeof(HF));
  if (!result.coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_new\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < len; i++) {
       result.coeffs[i] = coeffs[i];
  }
  poly_normalize(&result);
  return result;
}

POLY poly_zero() {
  HF zero = hf_new(0);
  return poly_new(&zero, 1);
}

POLY poly_one() {
  HF one = hf_new(1);
  return poly_new(&one, 1);
}

bool poly_is_zero(const POLY *polynomial) {
  return (polynomial->len == 1 && hf_equal(polynomial->coeffs[0], hf_new(0)));
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

void poly_divide(const POLY *numerator, const POLY *denominator, POLY *quotient, POLY *remainder) {
  // ensure the denominator is not zero
  if (denominator->len == 0 || hf_equal(denominator->coeffs[denominator->len - 1], hf_new(0))) {
    fprintf(stderr, "Division by zero polynomial in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // initialize quotient and remainder
  size_t num_len = numerator->len;
  size_t den_len = denominator->len;
  HF *quot_coeffs = (HF *)calloc(num_len, sizeof(HF));
  HF *rem_coeffs = (HF *)calloc(num_len, sizeof(HF));
  if (!quot_coeffs || !rem_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // copy numerator coefficients to remainder
  for (size_t i = 0; i < num_len; i++) {
    rem_coeffs[i] = numerator->coeffs[i];
  }

  HF den_lead_inv = hf_inv(denominator->coeffs[den_len - 1]);

  for (ssize_t i = num_len - 1; i >= (ssize_t)(den_len - 1); i--) {
    // compute the coefficient for the current term of the quotient
    HF coeff = hf_mul(rem_coeffs[i], den_lead_inv);
    quot_coeffs[i - (den_len - 1)] = coeff;

    // subtract the scaled denominator from the remainder
    for (ssize_t j = 0; j < (ssize_t)den_len; j++) {
      rem_coeffs[i - j] = hf_sub(rem_coeffs[i - j], hf_mul(coeff, denominator->coeffs[den_len - 1 - j]));
    }
  }

  // remove leading zeros from quotient and remainder
  size_t quot_len = num_len - den_len + 1;
  while (quot_len > 0 && hf_equal(quot_coeffs[quot_len - 1], hf_new(0))) {
    quot_len--;
  }

  size_t rem_len = den_len - 1;
  while (rem_len > 0 && hf_equal(rem_coeffs[rem_len - 1], hf_new(0))) {
    rem_len--;
  }

  // set the quotient and remainder polynomials
  *quotient = poly_new(quot_coeffs, quot_len);
  *remainder = poly_new(rem_coeffs, rem_len);

  // clean up temporary arrays
  free(quot_coeffs);
  free(rem_coeffs);
}


POLY poly_scale(const POLY *polynomial, HF scalar) {
  HF *coeffs = (HF *)calloc(polynomial->len, sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_scale\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < polynomial->len; i++) {
    coeffs[i] = hf_mul(polynomial->coeffs[i], scalar);
  }
  POLY result = poly_new(coeffs, polynomial->len);
  free(coeffs);
  return result;
}

POLY poly_shift(const POLY *polynomial, size_t shift) {
  // Multiply polynomial by x^shift
  size_t new_len = polynomial->len + shift;
  HF *new_coeffs = (HF *)calloc(new_len, sizeof(HF));
  if (!new_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_shift\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < polynomial->len; i++) {
    new_coeffs[i + shift] = polynomial->coeffs[i];
  }
  POLY result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

POLY poly_slice(const POLY *polynomial, size_t start, size_t end) {
  if (start >= end || end > polynomial->len) {
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
    new_coeffs[i] = polynomial->coeffs[start + i];
  }
  POLY result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

POLY poly_negate(const POLY *polynomial) {
  HF *neg_coeffs = (HF *)malloc(polynomial->len * sizeof(HF));
  if (!neg_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_negate\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < polynomial->len; i++) {
    neg_coeffs[i] = hf_neg(polynomial->coeffs[i]);
  }
  POLY result = poly_new(neg_coeffs, polynomial->len);
  free(neg_coeffs);
  return result;
}

void poly_free(POLY *polynomial) {
  if (polynomial->coeffs) {
    free(polynomial->coeffs);
    polynomial->coeffs = NULL;
  }
  polynomial->len = 0;
}

HF poly_eval(const POLY *polynomial, HF x) {
  HF x_pow = hf_one();
  HF y = hf_zero();
  for (size_t i = 0; i < polynomial->len; i++) {
    HF term = hf_mul(x_pow, polynomial->coeffs[i]);
    y = hf_add(y, term);
    x_pow = hf_mul(x_pow, x);
  }
  return y;
}

POLY poly_z(const HF *points, size_t len) {
  POLY acc = poly_one();
  for (size_t i = 0; i < len; i++) {
    HF neg_x = hf_neg(points[i]);
    HF coeffs[] = {neg_x, hf_new(1)};
    POLY term = poly_new(coeffs, 2);
    POLY temp = poly_mul(&acc, &term);
    poly_free(&acc);
    acc = temp;
  }
  return acc;
}

POLY poly_lagrange(const HF *x_points, const HF *y_points, size_t len) {
  POLY l = poly_zero();
  for (size_t j = 0; j < len; j++) {
    POLY l_j = poly_one();
    for (size_t i = 0; i < len; i++) {
      if (i != j) {
        HF denom = hf_sub(x_points[j], x_points[i]);
	HF denom_inv = hf_inv(denom);
	if (denom_inv.value == 0) {
	  fprintf(stderr, "Error: Lagrange polynomial x points must be unique\n");
	  exit(EXIT_FAILURE);
	}
	HF c = denom_inv;
	HF neg_cxi = hf_neg(hf_mul(c, x_points[i]));
	HF coeffs[] = {neg_cxi, c};
	POLY term = poly_new(coeffs, 2);
	POLY temp = poly_mul(&l_j, &term);
	poly_free(&l_j);
	l_j = temp;
	poly_free(&term);
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
