#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

#include "fe.h"
#include "poly.h"

typedef struct {
  size_t m; // number of rows
  size_t n; // number of cols
  u8_fe *v; // flat array of field elements
} matrix;

matrix matrix_zero(size_t m, size_t n) {
  matrix mat;
  mat.m = m;
  mat.n = n;
  mat.v = (u8_fe *)calloc(m * n, sizeof(u8_fe));
  if (!mat.v) {
    fprintf(stderr, "Memory allocation failed in matrix_zero\n");
    exit(EXIT_FAILURE);
  }
  return mat;
}

matrix matrix_new(u8_fe *v, size_t m, size_t n) {
  matrix mat;
  mat.m = m;
  mat.n = n;
  size_t size = m * n;
  mat.v = (u8_fe *)malloc(size * sizeof(u8_fe));
  if (!mat.v) {
    fprintf(stderr, "Memory allocation failed in matrix_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < size; i++)
    mat.v[i] = v[i];
  return mat;
}

// get element at (row, col)
u8_fe matrix_get(const matrix *mat, size_t row, size_t col) {
  if (row >= mat->m || col >= mat->n) {
    fprintf(stderr, "Index out of bounds in matrix_get\n");
    exit(EXIT_FAILURE);
  }
  return mat->v[col + row * mat->n];
}

// set element at (row, col)
void matrix_set(matrix *mat, size_t row, size_t col, u8_fe value) {
  if (row >= mat->m || col >= mat->n) {
    fprintf(stderr, "Index out of bounds in matrix_set\n");
    exit(EXIT_FAILURE);
  }
  mat->v[col + row * mat->n] = value;
}

void matrix_free(matrix *mat) {
  if (mat->v) {
    free(mat->v);
    mat->v = NULL;
  }
  mat->m = 0;
  mat->n = 0;
}

matrix matrix_add(const matrix *a, const matrix *b) {
  if (a->m != b->m || a->n != b->n) {
    fprintf(stderr, "Matrix dimensions must match for additoin\n");
    exit(EXIT_FAILURE);
  }
  matrix result = matrix_zero(a->m, a->n);
  size_t size = a->m * a->n;
  for (size_t i = 0; i < size; i++)
    result.v[i] = u8_fe_add(a->v[i], b->v[i]);
  return result;
}

matrix matrix_mul(const matrix *a, const matrix *b) {
  if (a->n != b->m) {
       fprintf(stderr, "Matrix multiplication error: Dimensions (%zu x %zu) and (%zu x %zu) incompatible.\n", a->m, a->n, b->m, b->n);
    exit(EXIT_FAILURE);
  }
  matrix result = matrix_zero(a->m, b->n);
  for (size_t i = 0; i < a->m; i++) {
    for (size_t j = 0; j < b->n; j++) {
      u8_fe sum = u8_fe_new(0);
      for (size_t k = 0; k < a->n; k++) {
        u8_fe product = u8_fe_mul(matrix_get(a, i, k), matrix_get(b, k, j));
        sum = u8_fe_add(sum, product);
      }
      matrix_set(&result ,i, j, sum);
    }
  }
  return result;
}

matrix matrix_mul_poly(const matrix *mat, const poly *p) {
  // convert polynomial coefficients to a column matrix
  matrix vec = matrix_zero(p->len, 1);
  for (size_t i = 0; i < p->len; i++)
    matrix_set(&vec, i, 0, p->coeffs[i]);
  // multiply the matrix by the vector
  matrix result = matrix_mul(mat, &vec);
  matrix_free(&vec);
  return result;
}

void matrix_gauss_jordan(matrix *mat) {
  size_t row_count = mat->m;
  size_t col_count = mat->n;
  size_t lead = 0;

  for (size_t r = 0; r < row_count; r++) {
    if (col_count <= lead)
      return;
    size_t i = r;
    while (u8_fe_equal(matrix_get(mat, i , lead), u8_fe_new(0))) {
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
        u8_fe temp = matrix_get(mat, i, k);
        matrix_set(mat, i, k, matrix_get(mat, r, k));
        matrix_set(mat, r, k, temp);
      }
    }
    // normalize row r
    u8_fe div = matrix_get(mat, r, lead);
    if (!u8_fe_equal(div, u8_fe_new(0))) {
      for (size_t k = 0; k < col_count; k++) {
        u8_fe value = matrix_get(mat, r, k);
        u8_fe result = u8_fe_div(value, div);
        matrix_set(mat, r, k, result);
      }
    }
    // eliminate other rows
    for (size_t i = 0; i < row_count; i++) {
      if (i != r) {
        u8_fe mult = matrix_get(mat, i, lead);
        for (size_t k = 0; k < col_count; k++) {
          u8_fe value = matrix_get(mat, i, k);
          u8_fe subtrahend = u8_fe_mul(matrix_get(mat, r, k), mult);
          u8_fe result = u8_fe_sub(value, subtrahend);;
          matrix_set(mat, i, k, result);
        }
      }
    }
    lead++;
  }
}

matrix matrix_inv(const matrix *mat) {
  if (mat->m != mat->n) {
    fprintf(stderr, "Only square matrices can be inverted\n");
    exit(EXIT_FAILURE);
  }
  size_t n = mat->n;
  // create augmented matrix [mat | I]
  matrix aug = matrix_zero(n, 2 * n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      matrix_set(&aug, i, j, matrix_get(mat, i, j));
    }
    matrix_set(&aug, i, i + n, u8_fe_new(1)); // identity matrix
  }
  // perform gauss-jordan elimination
  matrix_gauss_jordan(&aug);
  // extract the inverse matrix from the augmented matrix
  matrix inv = matrix_zero(n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      matrix_set(&inv, i, j, matrix_get(&aug, i, j + n));
    }
  }
  matrix_free(&aug);
  return inv;
}

#endif
