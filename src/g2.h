#ifndef G2_H
#define G2_H

#include <stdint.h>
#include "gf.h"

typedef struct {
  GF x, y;
} G2;

G2 g2_new(uint64_t x, uint64_t y) {
  G2 point = {
    .x = f101(x),
    .y = f101(y),
  };
  return point;
}

G2 g2_generator() {
  return g2_new(36, 31);
}

uint64_t g2_embedding_degree() {
  return 2;
}

G2 g2_neg(G2 *p) {
  GF neg_y = gf_neg(p->y);
  return g2_new(p->x.value, neg_y.value);
}

G2 g2_add(const G2 *p, const G2 *q) {
  GF x, y;

  if (gf_equal(p->x, q->x) && gf_equal(p->y, q->y)) {
    GF two = f101(2);
    GF three = f101(3);
    GF x1 = p->x;
    GF y1 = p->y;
    GF x_sq = gf_mul(x1, x1);
    GF num = gf_mul(three, x_sq);
    GF div = gf_mul(two, y1);
    GF m = gf_div(num, div);
    GF m_sq= gf_mul(m, m);
    GF neg_two = gf_neg(two);
    GF neg_two_inv = gf_inv(neg_two);
    GF m_sq_neg_two = gf_mul(m_sq, neg_two_inv);
    x = gf_sub(m_sq_neg_two, gf_mul(two, x1));
    y = gf_sub(gf_mul(gf_mul(neg_two_inv, m), gf_sub(gf_mul(three, x1), m_sq_neg_two)), y1);
  } else {
    GF x1 = p->x;
    GF y1 = p->y;
    GF x2 = q->x;
    GF y2 = q->y;
    GF num = gf_sub(y2, y1);
    GF div = gf_sub(x2, x1);
    GF m = gf_div(num, div);
    GF two = f101(2);
    GF neg_two = gf_neg(two);
    GF m_sq = gf_mul(m, m);
    GF m_sq_neg_two = gf_mul(m_sq, neg_two);
    x = gf_sub(gf_sub(m_sq_neg_two, x1), x2);
    y = gf_sub(gf_mul(m, gf_sub(x1, x)), y1);
  }
  return g2_new(x.value, y.value);
}

G2 g2_mul(G2 base, uint64_t scalar) {
     int flag = 0;
     G2 result;
     while (scalar > 0) {
       if (scalar % 2 == 1) {
	 if (flag) {
	   result = g2_add(&result, &base);
	 } else {
	   result = base;
	   flag = 1;
	 }
       }
       scalar >>= 1;
       base = g2_add(&base, &base);
     }
     return result;
}

#endif
