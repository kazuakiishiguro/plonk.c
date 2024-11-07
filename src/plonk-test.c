#include <assert.h>
#include "plonk.h"

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

int main() {
  test_plonk_init();
}
