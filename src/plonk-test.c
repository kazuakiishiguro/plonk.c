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

int main() {
  test_plonk_init();
  test_interpolate_at_h();
}
