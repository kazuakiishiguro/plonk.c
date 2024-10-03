#include <assert.h>
#include "pairing.h"

void test_pairing_vectors() {
  g1_p gen1 = g1_p_generator();
  g2_p gen2 = g2_p_generator();
  g1_p p = gen1;
  g1_p r = g1_p_mul(&gen1, 4);
  g2_p q = g2_p_mul(gen2, 3);
  GF a = gf_new(5);

  g1_p p_mul_a = g1_p_mul(&p, a.value);
  gtp left = pairing(&p_mul_a, &q);
  g2_p q_mul_a = g2_p_mul(q, a.value);
  gtp right = pairing(&p, &q_mul_a);
  assert(gtp_equal(&left, &right));

  gtp left2 = pairing(&p_mul_a, &q);
  gtp p_q = pairing(&p, &q);
  gtp right2 = gtp_pow(&p_q, a.value);
  assert(gtp_equal(&left2, &right2));

  g1_p tmp = g1_p_add(&p, &r);
  gtp p_plus_r = pairing(&tmp, &q);
  gtp r_q = pairing(&r, &q);
  gtp p_q_mul_r_q = gtp_mul(&p_q, &r_q);
  assert(gtp_equal(&p_plus_r, &p_q_mul_r_q));
}

int main() {
  test_pairing_vectors();
  return 0;
}
