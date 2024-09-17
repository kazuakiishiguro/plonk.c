#include <stdbool.h>
#include <stdint.h>
#include "fe.h"

typedef struct {
  u8_fe x, y;
  bool infinite;
} g1_p;

g1_p g1_p_new(uint64_t x_val, uint64_t y_val) {
  g1_p point = {
    .x = f101(x_val),
    .y = f101(y_val),
    .infinite = false,
  };
  return point;
}

g1_p g1_p_generator() {
  return g1_p_new(1, 2);
}

bool g1_p_is_on_curve(const g1_p *point) {
  if (point->infinite) return true;

  // y^2 = x^3 + 3 mod MODULO
  return u8_fe_equal(u8_fe_pow(point->y, 2), u8_fe_add(u8_fe_pow(point->x, 3), f101(3)));
}

g1_p g1_p_identity() {
  return (g1_p){.infinite = true};
}

g1_p g1_p_double(const g1_p *a) {
  if (a->infinite || u8_fe_equal(a->y, u8_fe_new(0))) return g1_p_identity();

  u8_fe two = f101(2);
  u8_fe three = f101(3);
  u8_fe x_sq = u8_fe_mul(a->x, a->x);
  u8_fe s = u8_fe_mul(three, x_sq);
  u8_fe d = u8_fe_mul(two, a->y);
  u8_fe m = u8_fe_div(s, d);

  u8_fe m_sq = u8_fe_mul(m, m);
  u8_fe x_r = u8_fe_mul(two, a->x);
  x_r = u8_fe_sub(m_sq, x_r);
  u8_fe y_r = u8_fe_mul(three, a->x);
  y_r = u8_fe_sub(y_r, m_sq);
  y_r = u8_fe_mul(m, y_r);
  y_r = u8_fe_sub(y_r, a->y);

  return g1_p_new(x_r.value, y_r.value);
}

// https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#G1P_addition
g1_p g1_p_add(const g1_p *a, const g1_p *b) {
  if (a->infinite) return *b;
  if (b->infinite) return *a;

  if (u8_fe_equal(a->x, b->x)) {
    u8_fe sum_y = u8_fe_add(a->y, b->y);
    if (u8_fe_equal(sum_y, u8_fe_new(0))) {
      // Pints are inverse of each other
      return g1_p_identity();
    } else {
      // Points are the same, perform doubling
      return g1_p_double(a);
    }
  }

  u8_fe s = u8_fe_sub(b->y, a->y);
  u8_fe d = u8_fe_sub(b->x, a->x);
  u8_fe m = u8_fe_mul(s, u8_fe_inv(d));

  u8_fe m_sq = u8_fe_mul(m, m);
  u8_fe x_r = u8_fe_sub(u8_fe_sub(m_sq, a->x), b->x);
  u8_fe y_r = u8_fe_sub(u8_fe_mul(m, u8_fe_sub(a->x, x_r)), a->y);

  return g1_p_new(x_r.value, y_r.value);
}

g1_p g1_p_neg(g1_p *a) {
  if (a->infinite) return *a;
  u8_fe neg_y = u8_fe_neg(a->y);
  return g1_p_new(a->x.value, neg_y.value);
}

g1_p g1_p_mul(const g1_p *point, uint64_t scalar) {
  g1_p result = g1_p_identity();
  g1_p added = *point;

  while(scalar > 0) {
    if (scalar & 1) {
      result = g1_p_add(&result, &added);
    }
    added = g1_p_double(&added);
    scalar >>= 1;
  }
  return result;
}
