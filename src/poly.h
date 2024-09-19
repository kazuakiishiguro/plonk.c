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

void poly_free(poly *p) {
  if (p->coeffs) {
    free(p->coeffs);
    p->coeffs = NULL;
  }
  p->len = 0;
}

u8_fe poly_eval(const poly *p, u8_fe x) {
  u8_fe x_pow = u8_fe_new(1);
  u8_fe y = p->coeffs[0];
  for (size_t i = 1; i < p->len; i++) {
    x_pow = u8_fe_mul(x_pow, x);
    u8_fe term = u8_fe_mul(x_pow, p->coeffs[i]);
    y = u8_fe_add(y, term);
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
