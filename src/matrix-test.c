#include <stdbool.h>
#include "matrix.h"

bool matrix_equal(const MATRIX *a, const MATRIX *b) {
  if (a->m != b->m || a->n != b->n)
    return false;
  size_t size = a->n * a->n;
  for (size_t i = 0; i < size; i++) {
    if (!hf_equal(a->v[i], b->v[i])) {
      return false;
    }
  }
  return true;
}

void test_matrix_add() {
  HF values_a[] = {f17(1), f17(2)};
  HF values_b[] = {f17(3), f17(4)};
  MATRIX a = matrix_new(values_a, 2, 1);
  MATRIX b = matrix_new(values_b, 2, 1);
  MATRIX result = matrix_add(&a, &b);
  HF values_expected[] = {f17(4), f17(6)};
  MATRIX expected = matrix_new(values_expected, 2, 1);
  matrix_equal(&result, &expected);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&result);
  matrix_free(&expected);
}

void test_matrix_mul() {
  HF values_a[] = {f17(1), f17(2), f17(3), f17(4), f17(5), f17(6)};
  HF values_b[] = {f17(10), f17(11), f17(20), f17(21), f17(30), f17(31)};
  MATRIX a = matrix_new(values_a, 2, 3);
  MATRIX b = matrix_new(values_b, 3, 2);
  MATRIX result = matrix_mul(&a, &b);
  HF values_expected[] = {f17(140), f17(146), f17(320), f17(335)};
  MATRIX expected = matrix_new(values_expected, 2, 2);
  matrix_equal(&result, &expected);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&result);
  matrix_free(&expected);
}

void test_matrix_inv() {
  HF values[] = {f17(1), f17(2), f17(3), f17(4)};
  MATRIX mat = matrix_new(values, 2, 2);
  MATRIX inv = matrix_inv(&mat);
  MATRIX inv_inv = matrix_inv(&inv);
  matrix_equal(&inv, &inv_inv);

  matrix_free(&mat);
  matrix_free(&inv);
  matrix_free(&inv_inv);
}

int main() {
  test_matrix_add();
  test_matrix_mul();
  test_matrix_inv();
  return 0;
}
