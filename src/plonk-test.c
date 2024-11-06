#include <assert.h>
#include "plonk.h"

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
}

int main() {
  test_plonk_init();
}
