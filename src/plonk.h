#ifndef PLONK_H
#define PLONK_H

#include <stdio.h>
#include <stdlib.h>

#include "pairing.h"

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

#endif // PLONK_H
