#ifndef PAIRING_H
#define PAIRING_H

#include <stdint.h>
#include "g1.h"
#include "g2.h"
#include "gt.h"

int gtp_equal(const gtp *p, const gtp *q) {
  return u8_fe_equal(p->a, q->a) && u8_fe_equal(p->b, q->b);
}

typedef struct {
  u8_fe x;
  u8_fe y;
  u8_fe c;
} line_eq;

line_eq line(const g1_p *a, const g1_p *b) {
  u8_fe m = u8_fe_sub(b->x, a->x);
  u8_fe n = u8_fe_sub(b->y, a->y);

  u8_fe x = n;
  u8_fe y = u8_fe_neg(m);
  u8_fe c = u8_fe_sub(u8_fe_mul(m, a->y), u8_fe_mul(n, a->x));
  line_eq l_eq = { .x = x, .y = y, .c = c };

  return l_eq;
}

gtp pairing_f(uint64_t r, const g1_p *p, const g2_p *q) {
  if (r == 1) {
    return gtp_new(f101(1), f101(0));
  } else if (r % 2 == 1) {
    uint64_t r_minus_1 = r - 1;
    g1_p rp = g1_p_mul(p, r_minus_1);
    line_eq l_eq = line(&rp, p);

    gtp prev = pairing_f(r_minus_1, p, q);

    u8_fe x = u8_fe_add(u8_fe_mul(q->x, l_eq.x), l_eq.c);
    u8_fe y = u8_fe_mul(q->y, l_eq.y);

    gtp gt_term = gtp_new(x, y);

    return gtp_mul(&prev, &gt_term);
  } else {
    uint64_t r_div_2 = r / 2;
    g1_p rp = g1_p_mul(p, r_div_2);
    g1_p neg_rp = g1_p_neg(&rp);
    g1_p two_neg_rp = g1_p_mul(&neg_rp, 2);
    line_eq l_eq = line(&rp, &two_neg_rp);

    gtp prev = pairing_f(r_div_2, p, q);
    gtp prev_sq = gtp_mul(&prev, &prev);

    u8_fe x = u8_fe_add(u8_fe_mul(q->x, l_eq.x), l_eq.c);
    u8_fe y = u8_fe_mul(q->y, l_eq.y);

    gtp gt_term = gtp_new(x, y);

    return gtp_mul(&prev_sq, &gt_term);
  }
}

gtp pairing(const g1_p *g1, const g2_p *g2) {
  uint64_t p_order = MODULO;
  uint64_t r = g1_p_generator_subgroup_size().value;
  uint64_t k = g2_p_embedding_degree();

  // Compute p^k mod r (if necessary)
  // Since p^k may be large, we need to handle this carefully
  uint64_t p_pow_k = 1;
  for (uint64_t i = 0; i < k; i++) {
    p_pow_k *= p_order;
  }

  uint64_t exp = (p_pow_k - 1) / r;

  gtp f = pairing_f(r, g1, g2);

  return gtp_pow(&f, exp);
}

#endif
