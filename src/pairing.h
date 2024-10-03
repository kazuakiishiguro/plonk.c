#ifndef PAIRING_H
#define PAIRING_H

#include <stdint.h>
#include "g1.h"
#include "g2.h"
#include "gt.h"

int gtp_equal(const GTP *p, const GTP *q) {
  return gf_equal(p->a, q->a) && gf_equal(p->b, q->b);
}

typedef struct {
  GF x;
  GF y;
  GF c;
} LINE_EQ;

LINE_EQ line(const G1 *a, const G1 *b) {
  GF m = gf_sub(b->x, a->x);
  GF n = gf_sub(b->y, a->y);

  GF x = n;
  GF y = gf_neg(m);
  GF c = gf_sub(gf_mul(m, a->y), gf_mul(n, a->x));
  LINE_EQ l_eq = { .x = x, .y = y, .c = c };

  return l_eq;
}

GTP pairing_f(uint64_t r, const G1 *p, const G2 *q) {
  if (r == 1) {
    return gtp_new(f101(1), f101(0));
  } else if (r % 2 == 1) {
    uint64_t r_minus_1 = r - 1;
    G1 rp = g1_mul(p, r_minus_1);
    LINE_EQ l_eq = line(&rp, p);

    GTP prev = pairing_f(r_minus_1, p, q);

    GF x = gf_add(gf_mul(q->x, l_eq.x), l_eq.c);
    GF y = gf_mul(q->y, l_eq.y);

    GTP gt_term = gtp_new(x, y);

    return gtp_mul(&prev, &gt_term);
  } else {
    uint64_t r_div_2 = r / 2;
    G1 rp = g1_mul(p, r_div_2);
    G1 neg_rp = g1_neg(&rp);
    G1 two_neg_rp = g1_mul(&neg_rp, 2);
    LINE_EQ l_eq = line(&rp, &two_neg_rp);

    GTP prev = pairing_f(r_div_2, p, q);
    GTP prev_sq = gtp_mul(&prev, &prev);

    GF x = gf_add(gf_mul(q->x, l_eq.x), l_eq.c);
    GF y = gf_mul(q->y, l_eq.y);

    GTP gt_term = gtp_new(x, y);

    return gtp_mul(&prev_sq, &gt_term);
  }
}

GTP pairing(const G1 *g1, const G2 *g2) {
  uint64_t p_order = MODULO;
  uint64_t r = g1_generator_subgroup_size().value;
  uint64_t k = g2_embedding_degree();

  // Compute p^k mod r (if necessary)
  // Since p^k may be large, we need to handle this carefully
  uint64_t p_pow_k = 1;
  for (uint64_t i = 0; i < k; i++) {
    p_pow_k *= p_order;
  }

  uint64_t exp = (p_pow_k - 1) / r;

  GTP f = pairing_f(r, g1, g2);

  return gtp_pow(&f, exp);
}

#endif
