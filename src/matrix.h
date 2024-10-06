#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

#include "hf.h"

typedef struct {
  size_t m; // number of rows
  size_t n; // number of cols
  HF *v; // flat array of field elements
} MATRIX;

MATRIX matrix_zero(size_t m, size_t n) {
  MATRIX result;
  result.m = m;
  result.n = n;
  result.v = (HF *)calloc(m * n, sizeof(HF));
  if (!result.v) {
    fprintf(stderr, "Memory allocation failed in matrix_zero\n");
    exit(EXIT_FAILURE);
  }
  return result;
}

MATRIX matrix_new(HF *v, size_t m, size_t n) {
  MATRIX result;
  result.m = m;
  result.n = n;
  size_t size = m * n;
  result.v = (HF *)malloc(size * sizeof(HF));
  if (!result.v) {
    fprintf(stderr, "Memory allocation failed in matrix_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < size; i++)
    result.v[i] = v[i];
  return result;
}

// get element at (row, col)
HF matrix_get(const MATRIX *matrix, size_t row, size_t col) {
  if (row >= matrix->m || col >= matrix->n) {
    fprintf(stderr, "Index out of bounds in matrix_get\n");
    exit(EXIT_FAILURE);
  }
  return matrix->v[col + row * matrix->n];
}

// set element at (row, col)
void matrix_set(MATRIX *matrix, size_t row, size_t col, HF value) {
  if (row >= matrix->m || col >= matrix->n) {
    fprintf(stderr, "Index out of bounds in matrix_set\n");
    exit(EXIT_FAILURE);
  }
  matrix->v[col + row * matrix->n] = value;
}

void matrix_free(MATRIX *matrix) {
  if (matrix->v) {
    free(matrix->v);
    matrix->v = NULL;
  }
  matrix->m = 0;
  matrix->n = 0;
}

MATRIX matrix_add(const MATRIX *a, const MATRIX *b) {
  if (a->m != b->m || a->n != b->n) {
    fprintf(stderr, "Matrix dimensions must match for additoin\n");
    exit(EXIT_FAILURE);
  }
  MATRIX result = matrix_zero(a->m, a->n);
  size_t size = a->m * a->n;
  for (size_t i = 0; i < size; i++)
    result.v[i] = hf_add(a->v[i], b->v[i]);
  return result;
}

MATRIX matrix_mul(const MATRIX *a, const MATRIX *b) {
  if (a->n != b->m) {
       fprintf(stderr, "Matrix multiplication error: Dimensions (%zu x %zu) and (%zu x %zu) incompatible.\n", a->m, a->n, b->m, b->n);
    exit(EXIT_FAILURE);
  }
  MATRIX result = matrix_zero(a->m, b->n);
  for (size_t i = 0; i < a->m; i++) {
    for (size_t j = 0; j < b->n; j++) {
      HF sum = hf_new(0);
      for (size_t k = 0; k < a->n; k++) {
        HF product = hf_mul(matrix_get(a, i, k), matrix_get(b, k, j));
        sum = hf_add(sum, product);
      }
      matrix_set(&result ,i, j, sum);
    }
  }
  return result;
}

void matrix_gauss_jordan(MATRIX *matrix) {
  size_t row_count = matrix->m;
  size_t col_count = matrix->n;
  size_t lead = 0;

  for (size_t r = 0; r < row_count; r++) {
    if (col_count <= lead)
      return;
    size_t i = r;
    while (hf_equal(matrix_get(matrix, i , lead), hf_new(0))) {
      i++;
      if (i == row_count) {
        i = r;
        lead++;
        if (lead == col_count)
          return;
      }
    }
    // swap rows i and r
    if (i != r) {
      for (size_t k = 0; k < col_count; k++) {
        HF temp = matrix_get(matrix, i, k);
        matrix_set(matrix, i, k, matrix_get(matrix, r, k));
        matrix_set(matrix, r, k, temp);
      }
    }
    // normalize row r
    HF div = matrix_get(matrix, r, lead);
    if (!hf_equal(div, hf_new(0))) {
      for (size_t k = 0; k < col_count; k++) {
        HF value = matrix_get(matrix, r, k);
        HF result = hf_div(value, div);
        matrix_set(matrix, r, k, result);
      }
    }
    // eliminate other rows
    for (size_t i = 0; i < row_count; i++) {
      if (i != r) {
        HF mult = matrix_get(matrix, i, lead);
        for (size_t k = 0; k < col_count; k++) {
          HF value = matrix_get(matrix, i, k);
          HF subtrahend = hf_mul(matrix_get(matrix, r, k), mult);
          HF result = hf_sub(value, subtrahend);;
          matrix_set(matrix, i, k, result);
        }
      }
    }
    lead++;
  }
}

MATRIX matrix_inv(const MATRIX *matrix) {
  if (matrix->m != matrix->n) {
    fprintf(stderr, "Only square matrices can be inverted\n");
    exit(EXIT_FAILURE);
  }
  size_t n = matrix->n;
  // create augmented MATRIX [mat | I]
  MATRIX aug = matrix_zero(n, 2 * n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      matrix_set(&aug, i, j, matrix_get(matrix, i, j));
    }
    matrix_set(&aug, i, i + n, hf_new(1)); // identity matrix
  }
  // perform gauss-jordan elimination
  matrix_gauss_jordan(&aug);
  // extract the inverse MATRIX from the augmented matrix
  MATRIX inv = matrix_zero(n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      matrix_set(&inv, i, j, matrix_get(&aug, i, j + n));
    }
  }
  matrix_free(&aug);
  return inv;
}

#endif
