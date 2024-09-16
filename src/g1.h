#include <stdbool.h>
#include <stdint.h>
#include "fe.h"

typedef struct {
  u64_fe x, y;
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
  u64_fe x_cube = u64_fe_pow(&point->x, 3);
  u64_fe b = f101(3);
  return u64_fe_pow(&point->y, 2).value == u64_fe_add(&x_cube, &b).value;
}

g1_p g1_p_identity() {
  return (g1_p){.infinite = true};
}

g1_p g1_p_double(const g1_p *a) {
  u64_fe two = f101(2);
  u64_fe three = f101(3);
  u64_fe x_sq = u64_fe_mul(&a->x, &a->x);
  u64_fe s = u64_fe_mul(&three, &x_sq);
  u64_fe d = u64_fe_mul(&two, &a->y);
  u64_fe m = u64_fe_div(&s, &d);

  u64_fe m_sq = u64_fe_mul(&m, &m);
  u64_fe x_r = u64_fe_mul(&two, &a->x);
  x_r = u64_fe_sub(&m_sq, &x_r);
  u64_fe y_r = u64_fe_mul(&three, &a->x);
  y_r = u64_fe_sub(&y_r, &m_sq);
  y_r = u64_fe_mul(&m, &y_r);
  y_r = u64_fe_sub(&y_r, &a->y);

  return g1_p_new(x_r.value, y_r.value);
}

g1_p g1_p_add_diff(const g1_p *a, const g1_p *b) {
  u64_fe s = u64_fe_sub(&b->y, &a->y);
  u64_fe d = u64_fe_sub(&b->x, &a->x);
  u64_fe m = u64_fe_div(&s, &d);

  u64_fe m_sq = u64_fe_mul(&m, &m);
  u64_fe x_r = u64_fe_sub(&m_sq, &a->x);
  x_r = u64_fe_sub(&x_r, &b->x);
  u64_fe y_r = u64_fe_sub(&a->x, &x_r);
  y_r = u64_fe_mul(&m, &y_r);
  y_r = u64_fe_sub(&y_r, &a->y);
  return g1_p_new(x_r.value, y_r.value);
}

// https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#G1P_addition
g1_p g1_p_add(const g1_p *a, const g1_p *b) {
  if (a->infinite) return *b;
  if (b->infinite) return *a;

  g1_p neg_b = g1_p_new(b->x.value, u64_fe_neg(&b->y).value);
  if (a->x.value == neg_b.x.value && a->y.value == neg_b.y.value) return g1_p_identity();

  if (a->x.value == b->x.value && a->y.value == b->y.value) {
    return g1_p_double(a);
  } else {
    return g1_p_add_diff(a, b);
  }
}

g1_p g1_p_neg(g1_p *a) {
  if (a->infinite) return *a;
  u64_fe neg_y = u64_fe_neg(&a->y);
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
