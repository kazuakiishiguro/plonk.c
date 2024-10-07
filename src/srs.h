#ifndef SRS_H
#define SRS_H

#include <stdio.h>
#include <stdlib.h>

#include "g1.h"
#include "g2.h"
#include "poly.h"

typedef struct {
  G1 *g1s;  // array of g1 points [1, s, s^2, ..., s^n]
  size_t len; // length of the g1s array
  G2 g2_1;  // g2 generator
  G2 g2_s;  // g2 generator * s
} SRS;

SRS srs_create(GF secret, size_t n) {
  SRS srs;
  srs.len = n + 1; // including the 0-th power
  srs.g1s = (G1 *)malloc(srs.len * sizeof(G1));
  if (!srs.g1s) {
    fprintf(stderr, "Mamory allocation failed in srs_create\n");
    exit(EXIT_FAILURE);
  }

  G1 id = g1_identity();

  // initialize g1s[0] = g1 generator
  srs.g1s[0] = id;
  // compute s^i for i = 1 to n
  GF s_pow = secret;
  for (size_t i = 0 ; i < srs.len; i++) {
    srs.g1s[i] = g1_mul(&id, s_pow.value);
    s_pow = gf_mul(s_pow, secret);
  }

  // initialize g2s
  srs.g2_1 = g2_generator();
  srs.g2_s = g2_mul(srs.g2_1, secret.value);

  return srs;
}

void srs_free(SRS *srs) {
  if (srs->g1s) {
    free(srs->g1s);
    srs->g1s = NULL;
  }
  srs->len = 0;
}

G1 srs_eval_at_s(const SRS *srs, const POLY *vs) {
  if (vs->len > srs->len) {
       fprintf(stderr, "Poynomial degree exceeds SRS size: POLY degree: %zu, SRS supports up to degree: %zu \n", vs->len, vs->len);
    exit(EXIT_FAILURE);
  }
  // compute a(s) = sum_{i=0}^{n-1} vs.coeffs[i] * g1s[i]
  G1 acc = g1_identity();
  for (size_t i = 0; i < vs->len; i++) {
    // multiply gls[i] by vs.coeffs[i]
    HF coeff = vs->coeffs[i];
    G1 term = g1_mul(&srs->g1s[i], coeff.value);
    acc = g1_add(&acc, &term);
  }

  return acc;
}

#endif // SRS_H
