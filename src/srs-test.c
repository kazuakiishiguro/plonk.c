#include <assert.h>
#include "srs.h"

#define ASSERT_G1(p1, p2)			      \
    do {					      \
        assert((p1).x.value == (p2).x.value);	      \
        assert((p1).y.value == (p2).y.value);	      \
    } while (0)

void test_srs() {
  GF secret = f101(5);
  size_t n = 5;
  SRS srs = srs_create(secret, n);
  assert(srs.len == n + 1);
  G1 g1_gen = g1_identity();

  ASSERT_G1(srs.g1s[0], g1_gen);

  HF poly_coeffs[] = {f17(1), f17(2), f17(3)};
  POLY poly = poly_new(poly_coeffs, 3);
  G1 eval = srs_eval_at_s(&srs, &poly);

  // expected eval
  // eval = g1s[0]*1 + g1s[1]*2 + g1s[2]*3
  G1 g1_0_1 = g1_mul(&srs.g1s[0], f101(1).value);
  G1 g1_1_2 = g1_mul(&srs.g1s[1], f101(2).value);
  G1 g1_2_3 = g1_mul(&srs.g1s[2], f101(3).value);
  G1 acc = g1_add(&g1_0_1, &g1_1_2);
  G1 expected = g1_add(&acc, &g1_2_3);

  ASSERT_G1(eval, expected);

  srs_free(&srs);
  poly_free(&poly);
}

int main() {
  test_srs();
  return 0;
}
