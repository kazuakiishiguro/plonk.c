#ifndef PLONK_H
#define PLONK_H

#include <stdlib.h>

#include "poly.h"
#include "srs.h"
#include "hf.h"

#define OMEGA_VALUE 4 // example value; should be a primitive roor of unity in the field

typedef struct {
  SRS srs;
  HF *h;        // array of HF elements
  size_t h_len;
  HF *k1_h;     // array of k1 * h elements
  HF *k2_h;     // array of k2 * h elements
  POLY z_h_x;   // polynomial z_h_x
} PLONK;

PLONK plonk_new(SRS srs, size_t n) {
  PLONK plonk;
  plonk.srs = srs;
  plonk.h_len = n;

  // init
  HF OMEGA = hf_new(OMEGA_VALUE);

  // compute h = [omega^0, omega^1, ...., omega^n]
  plonk.h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.h) {
    fprintf(stderr, "Memory allocation failed in plnk_new\n");
    exit(EXIT_FAILURE);
  }

  for (uint8_t i = 0; i < plonk.h_len; i++) {
    plonk.h[i] = hf_pow(OMEGA, i);
  }

  return plonk;
}

#endif // PLONK_H
