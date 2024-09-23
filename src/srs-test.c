#include <assert.h>
#include "srs.h"

#define ASSERT_G1(p1, p2)			      \
    do {					      \
        assert((p1).x.value == (p2).x.value);	      \
        assert((p1).y.value == (p2).y.value);	      \
    } while (0)

void test_srs() {
  u8_fe secret = f101(5);
  size_t n = 5;
  srs s = srs_create(secret, n);
  assert(s.len == n + 1);
  g1_p g1_gen = g1_p_identity();

  ASSERT_G1(s.g1s[0], g1_gen);

  u8_fe poly_coeffs[] = {f101(1), f101(2), f101(3)};
  poly p = poly_new(poly_coeffs, 3);
  g1_p eval = srs_eval_at_s(&s, &p);

  // expected eval
  // eval = g1s[0]*1 + g1s[1]*2 + g1s[2]*3
  g1_p g1_0_1 = g1_p_mul(&s.g1s[0], f101(1).value);
  g1_p g1_1_2 = g1_p_mul(&s.g1s[1], f101(2).value);
  g1_p g1_2_3 = g1_p_mul(&s.g1s[2], f101(3).value);
  g1_p acc = g1_p_add(&g1_0_1, &g1_1_2);
  g1_p expected = g1_p_add(&acc, &g1_2_3);

  ASSERT_G1(eval, expected);

  srs_free(&s);
  poly_free(&p);
}

int main() {
  test_srs();
  return 0;
}
