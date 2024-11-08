#include <assert.h>
#include "plonk.h"

#define ASSERT_POLY(p1, p2)                              \
  do {                                                   \
    assert(p1.len == p2.len);                            \
    for (size_t i = 0; i < p1.len; i++)                  \
      assert(hf_equal(p1.coeffs[i], p2.coeffs[i]));      \
  } while (0)

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

void test_plonk_init() {
  // create the trusted setup
  GF secret = f101(2); // the toxic waste
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);

  assert(plonk.h_len == h_len);

  // check roots of unity(H) generation
  for (uint8_t i = 0; i < h_len; i ++) {
    assert(hf_equal(plonk.h[i], hf_pow(hf_new(OMEGA_VALUE), i)));
  }

  HF vals[] = {f17(13), f17(13), f17(13), f17(13), f17(13), f17(16), f17(4), f17(1), f17(13), f17(4), f17(13), f17(4), f17(13), f17(1), f17(4), f17(16)};
  MATRIX expected = matrix_new(vals, h_len, h_len);
  assert(matrix_equal(&plonk.h_pows_inv, &expected));

  matrix_free(&expected);
  plonk_free(&plonk);
}

void test_interpolate_at_h() {
  GF secret = f101(2);
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);
  HF vec[] = {f17(3), f17(4), f17(0), f17(0)}; // adds paddings
  POLY result = interpolate_at_h(&plonk, vec, 4);
  // 6+x+4x^2+9x^3
  HF coeffs[] = {f17(6), f17(1), f17(4), f17(9)};
  POLY expected = poly_new(coeffs, 4);
  ASSERT_POLY(result, expected);
}

void test_copy_constraints() {
  GF secret = f101(2);
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);

  CONSTRAINTS constraints;
  constraints.num_constraints = h_len;
  constraints.c_a = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_b = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_c = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  if (!constraints.c_a || !constraints.c_b || !constraints.c_c) {
    fprintf(stderr, "Memory allocation fialed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // define copy constraints
  // c_a: [b(1), b(2), b(3), c(1)]
  constraints.c_a[0] = (COPY_OF){COPYOF_B, 1};
  constraints.c_a[1] = (COPY_OF){COPYOF_B, 2};
  constraints.c_a[2] = (COPY_OF){COPYOF_B, 3};
  constraints.c_a[3] = (COPY_OF){COPYOF_C, 1};

  // c_b: [a(1), a(2), a(3), c(2)]
  constraints.c_b[0] = (COPY_OF){COPYOF_A, 1};
  constraints.c_b[1] = (COPY_OF){COPYOF_A, 2};
  constraints.c_b[2] = (COPY_OF){COPYOF_A, 3};
  constraints.c_b[3] = (COPY_OF){COPYOF_C, 2};

  // c_c: [a(4), a(4), a(4), c(3)]
  constraints.c_c[0] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[1] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[2] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[3] = (COPY_OF){COPYOF_C, 3};

  HF *sigma_c_a = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  HF *sigma_c_b = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  HF *sigma_c_c = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  copy_constraints_to_roots(&plonk, constraints.c_a, constraints.num_constraints, sigma_c_a);
  copy_constraints_to_roots(&plonk, constraints.c_b, constraints.num_constraints, sigma_c_b);
  copy_constraints_to_roots(&plonk, constraints.c_c, constraints.num_constraints, sigma_c_c);

  HF expected_c_a[] = {f17(2), f17(8), f17(15), f17(3)};
  HF expected_c_b[] = {f17(1), f17(4), f17(16), f17(12)};
  HF expected_c_c[] = {f17(13), f17(9), f17(5), f17(14)};
  assert(hf_equal(*sigma_c_a, *expected_c_a));
  assert(hf_equal(*sigma_c_b, *expected_c_b));
  assert(hf_equal(*sigma_c_c, *expected_c_c));
}

int main() {
  test_plonk_init();
  test_interpolate_at_h();
  test_copy_constraints();
}
