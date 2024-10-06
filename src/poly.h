#ifndef POLY_H
#define POLY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include "gf.h"

typedef struct {
  GF *coeffs;
  size_t len;
} poly;

void poly_normalize(poly *p) {
  while (p->len > 1 && gf_equal(p->coeffs[p->len - 1], gf_new(0))) {
    p->len--;
  }
  p->coeffs = realloc(p->coeffs, p->len * sizeof(GF));
};

poly poly_new(GF *coeffs, size_t len) {
  poly p;
  p.len = len;
  p.coeffs = (GF *)malloc(len * sizeof(GF));
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
  GF zero = gf_new(0);
  return poly_new(&zero, 1);
}

poly poly_one() {
  GF one = gf_new(1);
  return poly_new(&one, 1);
}

bool poly_is_zero(const poly *p) {
  return (p->len == 1 && gf_equal(p->coeffs[0], gf_new(0)));
}

poly poly_add(const poly *a, const poly *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  GF *coeffs = (GF *)malloc(max_len * sizeof(GF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_add\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    GF a_coeff = (i < a->len) ? a->coeffs[i] : gf_new(0);
    GF b_coeff = (i < b->len) ? b->coeffs[i] : gf_new(0);
    coeffs[i] = gf_add(a_coeff, b_coeff);
  }
  poly p = poly_new(coeffs, max_len);
  free(coeffs);
  return p;
}

poly poly_sub(const poly *a, const poly *b) {
  size_t max_len = (a->len > b->len) ? a->len : b->len;
  GF *coeffs = (GF *)malloc(max_len * sizeof(GF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_sub\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < max_len; i++) {
    GF a_coeff = (i < a->len) ? a->coeffs[i] : gf_new(0);
    GF b_coeff = (i < b->len) ? b->coeffs[i] : gf_new(0);
    coeffs[i] = gf_sub(a_coeff, b_coeff);
  }
  poly p = poly_new(coeffs, max_len);
  free(coeffs);
  return p;
}

poly poly_mul(const poly *a, const poly *b) {
  size_t result_len = a->len + b->len - 1;
  GF *coeffs = (GF *)calloc(result_len, sizeof(GF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_mul\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < a->len; i++) {
    for (size_t j = 0; j < b->len; j++) {
      GF product = gf_mul(a->coeffs[i], b->coeffs[j]);
      coeffs[i + j] = gf_add(coeffs[i + j], product);
    }
  }
  poly p = poly_new(coeffs, result_len);
  free(coeffs);
  return p;
}

void poly_divide(const poly *numerator, const poly *denominator, poly *quotient, poly *remainder) {
  // ensure the denominator is not zero
  if (denominator->len == 0 || gf_equal(denominator->coeffs[denominator->len - 1], gf_new(0))) {
    fprintf(stderr, "Division by zero polynomial in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // initialize quotient and remainder
  size_t num_len = numerator->len;
  size_t den_len = denominator->len;
  GF *quot_coeffs = (GF *)calloc(num_len, sizeof(GF));
  GF *rem_coeffs = (GF *)calloc(num_len, sizeof(GF));
  if (!quot_coeffs || !rem_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_divide\n");
    exit(EXIT_FAILURE);
  }

  // copy numerator coefficients to remainder
  for (size_t i = 0; i < num_len; i++) {
    rem_coeffs[i] = numerator->coeffs[i];
  }

  GF den_lead_inv = gf_inv(denominator->coeffs[den_len - 1]);

  for (ssize_t i = num_len - 1; i >= (ssize_t)(den_len - 1); i--) {
    // compute the coefficient for the current term of the quotient
    GF coeff = gf_mul(rem_coeffs[i], den_lead_inv);
    quot_coeffs[i - (den_len - 1)] = coeff;

    // subtract the scaled denominator from the remainder
    for (ssize_t j = 0; j < (ssize_t)den_len; j++) {
      rem_coeffs[i - j] = gf_sub(rem_coeffs[i - j], gf_mul(coeff, denominator->coeffs[den_len - 1 - j]));
    }
  }

  // remove leading zeros from quotient and remainder
  size_t quot_len = num_len - den_len + 1;
  while (quot_len > 0 && gf_equal(quot_coeffs[quot_len - 1], gf_new(0))) {
    quot_len--;
  }

  size_t rem_len = den_len - 1;
  while (rem_len > 0 && gf_equal(rem_coeffs[rem_len - 1], gf_new(0))) {
    rem_len--;
  }

  // set the quotient and remainder polynomials
  *quotient = poly_new(quot_coeffs, quot_len);
  *remainder = poly_new(rem_coeffs, rem_len);

  // clean up temporary arrays
  free(quot_coeffs);
  free(rem_coeffs);
}


poly poly_scale(const poly *p, GF scalar) {
  GF *coeffs = (GF *)calloc(p->len, sizeof(GF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_scale\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->len; i++) {
    coeffs[i] = gf_mul(p->coeffs[i], scalar);
  }
  poly result = poly_new(coeffs, p->len);
  free(coeffs);
  return result;
}

poly poly_shift(const poly *p, size_t shift) {
  // Multiply polynomial by x^shift
  size_t new_len = p->len + shift;
  GF *new_coeffs = (GF *)calloc(new_len, sizeof(GF));
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
  GF *new_coeffs = (GF *)malloc(new_len * sizeof(GF));
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
  GF *neg_coeffs = (GF *)malloc(p->len * sizeof(GF));
  if (!neg_coeffs) {
    fprintf(stderr, "Memory allocation failed in poly_negate\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->len; i++) {
    neg_coeffs[i] = gf_neg(p->coeffs[i]);
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

GF poly_eval(const poly *p, GF x) {
  GF x_pow = gf_one();
  GF y = gf_zero();
  for (size_t i = 0; i < p->len; i++) {
    GF term = gf_mul(x_pow, p->coeffs[i]);
    y = gf_add(y, term);
    x_pow = gf_mul(x_pow, x);
  }
  return y;
}

poly poly_z(const GF *points, size_t len) {
  poly acc = poly_one();
  for (size_t i = 0; i < len; i++) {
    GF neg_x = gf_neg(points[i]);
    GF coeffs[] = {neg_x, gf_new(1)};
    poly term = poly_new(coeffs, 2);
    poly temp = poly_mul(&acc, &term);
    poly_free(&acc);
    acc = temp;
  }
  return acc;
}

poly poly_lagrange(const GF *x_points, const GF *y_points, size_t len) {
  poly l = poly_zero();
  for (size_t j = 0; j < len; j++) {
    poly l_j = poly_one();
    for (size_t i = 0; i < len; i++) {
      if (i != j) {
        GF denom = gf_sub(x_points[j], x_points[i]);
	GF denom_inv = gf_inv(denom);
	if (denom_inv.value == 0) {
	  fprintf(stderr, "Error: Lagrange polynomial x points must be unique\n");
	  exit(EXIT_FAILURE);
	}
	GF c = denom_inv;
	GF neg_cxi = gf_neg(gf_mul(c, x_points[i]));
	GF coeffs[] = {neg_cxi, c};
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
