#ifndef PLONK_H
#define PLONK_H

#include <stdio.h>
#include <stdlib.h>

#include "pairing.h"
#include "poly.h"

typedef struct {
  // fields
  u8_fe gf;
  u8_fe hf;

  // ec points
  g1_p g1_generator;
  g2_p g2_generator;

  // constants
  u8_fe k1; // k1 x omega coset generator
  u8_fe k2; // k2 x omega coset generator
  u8_fe omega; // the generator in hf

  // function to map hf to gf
  u8_fe (*g_f)(u8_fe);
} plonk_types;

typedef struct {
  g1_p *g1s;
  g2_p g2_1; // g2 generator
  g2_p g2_s; // g2 generator * s
  size_t n; // number of g1 points
  plonk_types *pt;
} srs;

void srs_create(srs *s, plonk_types *pt, u8_fe sec, size_t n) {
  s->pt = pt;
  s->n = n;
  s->g1s = (g1_p *)malloc(n * sizeof(g1_p));
  if (!s->g1s) {
    fprintf(stderr, "Mamory allocation failed in srs_create\n");
    exit(EXIT_FAILURE);
  }

  // initialize g1s
  s->g1s[0] = g1_p_identity();
  for (size_t i = 0 ; i < n; i++) {
    s->g1s[i] = g1_p_mul(&s->g1s[i - 1], sec.value);
  }

  // initialize g2s
  s->g2_1 = g2_p_generator();
  s->g2_s = g2_p_mul(s->pt->g2_generator, sec.value);
}

void srs_clear(srs *s) {
  if (s->g1s) {
    /* for (size_t i = 0; i < s->n; i++) */
    /*   g1_p_clear(&s->g1s[i]); */
    free(s->g1s);
    s->g1s = NULL;
  }
}

g1_p srs_eval_at_s(srs *s, poly *vs) {
  // compute a(s) = sum_{i=0}^{n-1} vs.coeffs[i] * g1s[i]
  g1_p a_s = g1_p_identity();

  for (size_t i = 0; i < vs->len && i < s->n; i++) {
    // multiply gls[i] by vs.coeffs[i]
    g1_p term = g1_p_mul(&s->g1s[i], (uint64_t)&vs->coeffs[i].value);
    a_s = g1_p_add(&a_s, &term);
  }

  return a_s;
}

#endif // PLONK_H
