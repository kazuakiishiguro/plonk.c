#ifndef PLONK_H
#define PLONK_H

#include <stdlib.h>
#include "constraints.h"
#include "matrix.h"
#include "poly.h"
#include "srs.h"
#include "hf.h"

#define OMEGA_VALUE 4 // example value; should be a primitive roor of unity in the field
#define K1_VALUE    2 // should not be in H
#define K2_VALUE    3 // should not be in H or K1 * H

typedef struct {
  SRS srs;
  HF *h;        // array of HF elements
  MATRIX h_pows_inv;
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
  HF K1 = hf_new(K1_VALUE);
  HF K2 = hf_new(K2_VALUE);

  // compute h = [omega^0, omega^1, ...., omega^n]
  plonk.h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.h) {
    fprintf(stderr, "Memory allocation failed in plnk_new\n");
    exit(EXIT_FAILURE);
  }
  for (uint8_t i = 0; i < plonk.h_len; i++) {
    plonk.h[i] = hf_pow(OMEGA, i);
  }

  // check K1 and K2 are not in H
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.h[i], K1) || hf_equal(plonk.h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }

  // compute h1_h: h1_h[i] = K1 * h[i]
  // compute k2_h: k2_h[i] = K2 * h[i]
  plonk.k1_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  plonk.k2_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.k1_h || !plonk.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  // k1_h is a coset of H
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k1_h[i] = hf_mul(plonk.h[i], K1);
  }
  // check K2 is chosen so that it is neither an element of H nor k1_h[i]
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.k1_h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }
  // k2_h is a coset of H
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k2_h[i] = hf_mul(plonk.h[i], K2);
  }

  // build h_pows matrix (vandermonde matrix) and compute its inverse
  MATRIX h_pows = matrix_zero(plonk.h_len, plonk.h_len);
  for (size_t c = 0; c < h_pows.n; c++) {
    for (size_t r = 0; r < h_pows.m; r++) {
      matrix_set(&h_pows, r, c, hf_pow(plonk.h[r], c));
    }
  }
  plonk.h_pows_inv = matrix_inv(&h_pows);
  matrix_free(&h_pows);

  // compute z_h_x = Poly::z(h)
  plonk.z_h_x = poly_z(plonk.h, plonk.h_len);

  return plonk;
}

void plonk_free(PLONK *plonk) {
  srs_free(&plonk->srs);
  matrix_free(&plonk->h_pows_inv);
  poly_free(&plonk->z_h_x);

  if (plonk->h) {
    free(plonk->h);
    plonk->h = NULL;
  }

  if (plonk->k1_h) {
    free(plonk->k1_h);
    plonk->k1_h = NULL;
  }

  if (plonk->k2_h) {
    free(plonk->k2_h);
    plonk->k2_h = NULL;
  }
}

void copy_constraints_to_roots(const PLONK *plonk, const COPY_OF *copy_of, size_t len, HF *sigma) {
  for (size_t i = 0; i < len; i++) {
    size_t idx = copy_of[i].index - 1;
    switch (copy_of[i].type) {
    case COPYOF_A:
      sigma[i] = plonk->h[idx];
      break;
    case COPYOF_B:
      sigma[i] = plonk->k1_h[idx];
      break;
    case COPYOF_C:
      sigma[i] = plonk->k2_h[idx];
      break;
    default:
      fprintf(stderr, "Invalid copy_of type\n");
      exit(EXIT_FAILURE);
    }
  }
}

POLY interpolate_at_h(const PLONK *plonk, const HF *values, size_t len) {
  // ensure len == plonk->h_len
  if (len != plonk->h_len) {
    fprintf(stderr, "Length mismatch in interpolate_at_h: len= %zu, plonk->h_len = %zu\n", len, plonk->h_len);
    exit(EXIT_FAILURE);
  }

  MATRIX vec = matrix_zero(len, 1);
  for (size_t i = 0; i < len; i++) {
    matrix_set(&vec, i, 0, values[i]);
  }

  // multiply h_pows_inv * vec
  MATRIX res = matrix_mul(&plonk->h_pows_inv, &vec);

  // convert result matrix to polynomial
  HF *coeffs = (HF *)malloc(res.m * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }

  POLY result = poly_new(coeffs, res.m);

  // clean up
  matrix_free(&vec);
  matrix_free(&res);
  free(coeffs);

  return result;
}

#endif // PLONK_H
