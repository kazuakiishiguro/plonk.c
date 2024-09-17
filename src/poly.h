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

void poly_free(poly *p) {
  if (p->coeffs) {
    free(p->coeffs);
    p->coeffs = NULL;
  }
  p->len = 0;
}
