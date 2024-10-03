#include <assert.h>
#include "pairing.h"

void test_pairing_vectors() {
  G1 gen1 = g1_generator();
  G2 gen2 = g2_generator();
  G1 p = gen1;
  G1 r = g1_mul(&gen1, 4);
  G2 q = g2_mul(gen2, 3);
  GF a = gf_new(5);

  G1 p_mul_a = g1_mul(&p, a.value);
  GTP left = pairing(&p_mul_a, &q);
  G2 q_mul_a = g2_mul(q, a.value);
  GTP right = pairing(&p, &q_mul_a);
  assert(gtp_equal(&left, &right));

  GTP left2 = pairing(&p_mul_a, &q);
  GTP p_q = pairing(&p, &q);
  GTP right2 = gtp_pow(&p_q, a.value);
  assert(gtp_equal(&left2, &right2));

  G1 tmp = g1_add(&p, &r);
  GTP p_plus_r = pairing(&tmp, &q);
  GTP r_q = pairing(&r, &q);
  GTP p_q_mul_r_q = gtp_mul(&p_q, &r_q);
  assert(gtp_equal(&p_plus_r, &p_q_mul_r_q));
}

int main() {
  test_pairing_vectors();
  return 0;
}
