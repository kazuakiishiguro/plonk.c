#include <stdbool.h>
#include "matrix.h"

bool matrix_equal(const matrix *a, const matrix *b) {
  if (a->m != b->m || a->n != b->n)
    return false;
  size_t size = a->n * a->n;
  for (size_t i = 0; i < size; i++) {
    if (!u8_fe_equal(a->v[i], b->v[i])) {
      return false;
    }
  }
  return true;
}

void test_matrix_add() {
  u8_fe values_a[] = {f101(1), f101(2)};
  u8_fe values_b[] = {f101(3), f101(4)};
  matrix a = matrix_new(values_a, 2, 1);
  matrix b = matrix_new(values_b, 2, 1);
  matrix result = matrix_add(&a, &b);
  u8_fe values_expected[] = {f101(4), f101(6)};
  matrix expected = matrix_new(values_expected, 2, 1);
  matrix_equal(&result, &expected);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&result);
  matrix_free(&expected);
}

int main() {
  test_matrix_add();
  return 0;
}
