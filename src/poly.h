#ifndef POLY_H
#define POLY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include "fe.h"

typedef struct {
  u8_fe *coeffs;
  size_t len;
} poly;

void poly_normalize(poly *p) {
  while (p->len > 1 && u8_fe_equal(p->coeffs[p->len - 1], u8_fe_new(0))) {
    p->len--;
  }
  p->coeffs = realloc(p->coeffs, p->len * sizeof(u8_fe));
};

poly poly_new(u8_fe *coeffs, size_t len) {
  poly p;
  p.len = len;
  p.coeffs = (u8_fe *)malloc(len * sizeof(u8_fe));
  if (!p.coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_new\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < len; i++) {
       p.coeffs[i] = coeffs[i];
  }
  poly_normalize(&p);
  return p;
}

poly poly_zero() {
  u8_fe zero = u8_fe_new(0);
  return poly_new(&zero, 1);
}

poly poly_one() {
  u8_fe one = u8_fe_new(1);
  return poly_new(&one, 1);
}

bool poly_is_zero(const poly *p) {
  return (p->len == 1 && u8_fe_equal(p->coeffs[0], u8_fe_new(0)));
}

poly poly_add(const poly *a, const poly *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  u8_fe *coeffs = (u8_fe *)malloc(max_len * sizeof(u8_fe));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_add\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    u8_fe a_coeff = (i < a->len) ? a->coeffs[i] : u8_fe_new(0);
    u8_fe b_coeff = (i < b->len) ? b->coeffs[i] : u8_fe_new(0);
    coeffs[i] = u8_fe_add(a_coeff, b_coeff);
  }
  poly p = poly_new(coeffs, max_len);
  free(coeffs);
  return p;
}

poly poly_sub(const poly *a, const poly *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  u8_fe *coeffs = (u8_fe *)malloc(max_len * sizeof(u8_fe));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_sub\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    u8_fe a_coeff = (i < a->len) ? a->coeffs[i] : u8_fe_new(0);
    u8_fe b_coeff = (i < b->len) ? b->coeffs[i] : u8_fe_new(0);
    coeffs[i] = u8_fe_sub(a_coeff, b_coeff);
  }
  poly p = poly_new(coeffs, max_len);
  free(coeffs);
  return p;
}

poly poly_mul(const poly *a, const poly *b) {
  size_t result_len = a->len + b->len - 1;
  u8_fe *coeffs = (u8_fe *)calloc(result_len, sizeof(u8_fe));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_mul\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < a->len; i++) {
    for (size_t j = 0; j < b->len; j++) {
      u8_fe product = u8_fe_mul(a->coeffs[i], b->coeffs[j]);
      coeffs[i + j] = u8_fe_add(coeffs[i + j], product);
    }
  }
  poly p = poly_new(coeffs, result_len);
  free(coeffs);
  return p;
}

void poly_divide(const poly *numerator, const poly *denominator, poly *quotient, poly *remainder) {
  // ensure the denominator is not zero
  if (denominator->len == 0 || u8_fe_equal(denominator->coeffs[denominator->len - 1], u8_fe_new(0))) {
    fprintf(stderr, "Division by zero polynomial in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // initialize quotient and remainder
  size_t num_len = numerator->len;
  size_t den_len = denominator->len;
  u8_fe *quot_coeffs = (u8_fe *)calloc(num_len, sizeof(u8_fe));
  u8_fe *rem_coeffs = (u8_fe *)calloc(num_len, sizeof(u8_fe));
  if (!quot_coeffs || !rem_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // copy numerator coefficients to remainder
  for (size_t i = 0; i < num_len; i++) {
    rem_coeffs[i] = numerator->coeffs[i];
  }

  u8_fe den_lead_inv = u8_fe_inv(denominator->coeffs[den_len - 1]);

  for (ssize_t i = num_len - 1; i >= (ssize_t)(den_len - 1); i--) {
    // compute the coefficient for the current term of the quotient
    u8_fe coeff = u8_fe_mul(rem_coeffs[i], den_lead_inv);
    quot_coeffs[i - (den_len - 1)] = coeff;

    // subtract the scaled denominator from the remainder
    for (ssize_t j = 0; j < (ssize_t)den_len; j++) {
      rem_coeffs[i - j] = u8_fe_sub(rem_coeffs[i - j], u8_fe_mul(coeff, denominator->coeffs[den_len - 1 - j]));
    }
  }

  // remove leading zeros from quotient and remainder
  size_t quot_len = num_len - den_len + 1;
  while (quot_len > 0 && u8_fe_equal(quot_coeffs[quot_len - 1], u8_fe_new(0))) {
    quot_len--;
  }

  size_t rem_len = den_len - 1;
  while (rem_len > 0 && u8_fe_equal(rem_coeffs[rem_len - 1], u8_fe_new(0))) {
    rem_len--;
  }

  // set the quotient and remainder polynomials
  *quotient = poly_new(quot_coeffs, quot_len);
  *remainder = poly_new(rem_coeffs, rem_len);

  // clean up temporary arrays
  free(quot_coeffs);
  free(rem_coeffs);
}


poly poly_scale(const poly *p, u8_fe scalar) {
  u8_fe *coeffs = (u8_fe *)calloc(p->len, sizeof(u8_fe));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_scale\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->len; i++) {
    coeffs[i] = u8_fe_mul(p->coeffs[i], scalar);
  }
  poly result = poly_new(coeffs, p->len);
  free(coeffs);
  return result;
}

poly poly_shift(const poly *p, size_t shift) {
  // Multiply polynomial by x^shift
  size_t new_len = p->len + shift;
  u8_fe *new_coeffs = (u8_fe *)calloc(new_len, sizeof(u8_fe));
  if (!new_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_shift\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->len; i++) {
    new_coeffs[i + shift] = p->coeffs[i];
  }
  poly result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

poly poly_slice(const poly *p, size_t start, size_t end) {
  if (start >= end || end > p->len) {
    fprintf(stderr, "Invalid slice indices in poly_slice\n");
    exit(EXIT_FAILURE);
  }
  size_t new_len = end - start;
  u8_fe *new_coeffs = (u8_fe *)malloc(new_len * sizeof(u8_fe));
  if (!new_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_slice\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < new_len; i++) {
    new_coeffs[i] = p->coeffs[start + i];
  }
  poly result = poly_new(new_coeffs, new_len);
  free(new_coeffs);
  return result;
}

poly poly_negate(const poly *p) {
  u8_fe *neg_coeffs = (u8_fe *)malloc(p->len * sizeof(u8_fe));
  if (!neg_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_negate\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->len; i++) {
    neg_coeffs[i] = u8_fe_neg(p->coeffs[i]);
  }
  poly result = poly_new(neg_coeffs, p->len);
  free(neg_coeffs);
  return result;
}

void poly_free(poly *p) {
  if (p->coeffs) {
    free(p->coeffs);
    p->coeffs = NULL;
  }
  p->len = 0;
}

u8_fe poly_eval(const poly *p, hf_fe x_hf) {
  u8_fe x = gf_from_hf(x_hf);
  u8_fe x_pow = u8_fe_one();
  u8_fe y = u8_fe_zero();
  for (size_t i = 0; i < p->len; i++) {
    u8_fe term = u8_fe_mul(x_pow, p->coeffs[i]);
    y = u8_fe_add(y, term);
    x_pow = u8_fe_mul(x_pow, x);
  }
  return y;
}

poly poly_z(const u8_fe *points, size_t len) {
  poly acc = poly_one();
  for (size_t i = 0; i < len; i++) {
    u8_fe neg_x = u8_fe_neg(points[i]);
    u8_fe coeffs[] = {neg_x, u8_fe_new(1)};
    poly term = poly_new(coeffs, 2);
    poly temp = poly_mul(&acc, &term);
    poly_free(&acc);
    acc = temp;
  }
  return acc;
}

poly poly_lagrange(const u8_fe *x_points, const u8_fe *y_points, size_t len) {
  poly l = poly_zero();
  for (size_t j = 0; j < len; j++) {
    poly l_j = poly_one();
    for (size_t i = 0; i < len; i++) {
      if (i != j) {
        u8_fe denom = u8_fe_sub(x_points[j], x_points[i]);
	u8_fe denom_inv = u8_fe_inv(denom);
	if (denom_inv.value == 0) {
	  fprintf(stderr, "Error: Lagrange polynomial x points must be unique\n");
	  exit(EXIT_FAILURE);
	}
	u8_fe c = denom_inv;
	u8_fe neg_cxi = u8_fe_neg(u8_fe_mul(c, x_points[i]));
	u8_fe coeffs[] = {neg_cxi, c};
	poly term = poly_new(coeffs, 2);
	poly temp = poly_mul(&l_j, &term);
	poly_free(&l_j);
	l_j = temp;
	poly_free(&term);
      }
    }
    poly scaled_l_j = poly_scale(&l_j, y_points[j]);
    poly temp_l = poly_add(&l, &scaled_l_j);
    poly_free(&l);
    poly_free(&l_j);
    poly_free(&scaled_l_j);
    l = temp_l;
  }
  return l;
}

#endif // POLY_H
