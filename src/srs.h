#ifndef SRS_H
#define SRS_H

#include <stdio.h>
#include <stdlib.h>

#include "g1.h"
#include "g2.h"
#include "poly.h"

typedef struct {
  g1_p *g1s;  // array of g1 points [1, s, s^2, ..., s^n]
  size_t len; // length of the g1s array
  g2_p g2_1;  // g2 generator
  g2_p g2_s;  // g2 generator * s
} srs;

srs srs_create(u8_fe secret, size_t n) {
  srs s;
  s.len = n + 1; // including the 0-th power
  s.g1s = (g1_p *)malloc(s.len * sizeof(g1_p));
  if (!s.g1s) {
    fprintf(stderr, "Mamory allocation failed in srs_create\n");
    exit(EXIT_FAILURE);
  }

  g1_p id = g1_p_identity();

  // initialize g1s[0] = g1 generator
  s.g1s[0] = id;
  // compute s^i for i = 1 to n
  u8_fe s_pow = secret;
  for (size_t i = 0 ; i < s.len; i++) {
    s.g1s[i] = g1_p_mul(&id, s_pow.value);
    s_pow = u8_fe_mul(s_pow, secret);
  }

  // initialize g2s
  s.g2_1 = g2_p_generator();
  s.g2_s = g2_p_mul(s.g2_1, secret.value);

  return s;
}

void srs_free(srs *s) {
  if (s->g1s) {
    free(s->g1s);
    s->g1s = NULL;
  }
  s->len = 0;

}

g1_p srs_eval_at_s(const srs *s, const poly *vs) {
  // compute a(s) = sum_{i=0}^{n-1} vs.coeffs[i] * g1s[i]
  g1_p acc = g1_p_identity();
  size_t len = vs->len;
  if (len > s->len) {
       fprintf(stderr, "Poynomial degree exceeds SRS size: len: %zu, s->len: %zu \n", len, vs->len);
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < len; i++) {
    // multiply gls[i] by vs.coeffs[i]
    u8_fe coeff = vs->coeffs[i];
    g1_p term = g1_p_mul(&s->g1s[i], coeff.value);
    acc = g1_p_add(&acc, &term);
  }

  return acc;
}

#endif // SRS_H
