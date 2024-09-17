#ifndef G2_H
#define G2_H

#include <stdint.h>
#include "fe.h"

typedef struct {
  u8_fe x, y;
} g2_p;

g2_p g2_p_new(uint64_t x, uint64_t y) {
  g2_p point = {
    .x = f101(x),
    .y = f101(y),
  };
  return point;
}

g2_p g2_p_generator() {
  return g2_p_new(36, 31);
}

uint64_t g2_p_embedding_degree() {
  return 2;
}

g2_p g2_p_neg(g2_p *p) {
  u8_fe neg_y = u8_fe_neg(p->y);
  return g2_p_new(p->x.value, neg_y.value);
}

g2_p g2_p_add(const g2_p *p, const g2_p *q) {
  u8_fe x, y;

  if (u8_fe_equal(p->x, q->x) && u8_fe_equal(p->y, q->y)) {
    u8_fe two = f101(2);
    u8_fe three = f101(3);
    u8_fe x1 = p->x;
    u8_fe y1 = p->y;
    u8_fe x_sq = u8_fe_mul(x1, x1);
    u8_fe num = u8_fe_mul(three, x_sq);
    u8_fe div = u8_fe_mul(two, y1);
    u8_fe m = u8_fe_div(num, div);
    u8_fe m_sq= u8_fe_mul(m, m);
    u8_fe neg_two = u8_fe_neg(two);
    u8_fe neg_two_inv = u8_fe_inv(neg_two);
    u8_fe m_sq_neg_two = u8_fe_mul(m_sq, neg_two_inv);
    x = u8_fe_sub(m_sq_neg_two, u8_fe_mul(two, x1));
    y = u8_fe_sub(u8_fe_mul(u8_fe_mul(neg_two_inv, m), u8_fe_sub(u8_fe_mul(three, x1), m_sq_neg_two)), y1);
  } else {
    u8_fe x1 = p->x;
    u8_fe y1 = p->y;
    u8_fe x2 = q->x;
    u8_fe y2 = q->y;
    u8_fe num = u8_fe_sub(y2, y1);
    u8_fe div = u8_fe_sub(x2, x1);
    u8_fe m = u8_fe_div(num, div);
    u8_fe two = f101(2);
    u8_fe neg_two = u8_fe_neg(two);
    u8_fe m_sq = u8_fe_mul(m, m);
    u8_fe m_sq_neg_two = u8_fe_mul(m_sq, neg_two);
    x = u8_fe_sub(u8_fe_sub(m_sq_neg_two, x1), x2);
    y = u8_fe_sub(u8_fe_mul(m, u8_fe_sub(x1, x)), y1);
  }
  return g2_p_new(x.value, y.value);
}

g2_p g2_p_mul(g2_p base, uint64_t scalar) {
     int flag = 0;
     g2_p result;
     while (scalar > 0) {
       if (scalar % 2 == 1) {
	 if (flag) {
	   result = g2_p_add(&result, &base);
	 } else {
	   result = base;
	   flag = 1;
	 }
       }
       scalar >>= 1;
       base = g2_p_add(&base, &base);
     }
     return result;
}

#endif
