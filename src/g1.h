#ifndef G1_H
#define G1_H

#include <stdbool.h>
#include <stdint.h>
#include "gf.h"

typedef struct {
  GF x, y;
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
  return gf_equal(gf_pow(point->y, 2), gf_add(gf_pow(point->x, 3), f101(3)));
}

g1_p g1_p_identity() {
  return (g1_p){.infinite = true};
}

g1_p g1_p_double(const g1_p *a) {
  if (a->infinite || gf_equal(a->y, gf_new(0))) return g1_p_identity();

  GF two = f101(2);
  GF three = f101(3);
  GF x_sq = gf_mul(a->x, a->x);
  GF s = gf_mul(three, x_sq);
  GF d = gf_mul(two, a->y);
  GF m = gf_div(s, d);

  GF m_sq = gf_mul(m, m);
  GF x_r = gf_mul(two, a->x);
  x_r = gf_sub(m_sq, x_r);
  GF y_r = gf_mul(three, a->x);
  y_r = gf_sub(y_r, m_sq);
  y_r = gf_mul(m, y_r);
  y_r = gf_sub(y_r, a->y);

  return g1_p_new(x_r.value, y_r.value);
}

// https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#G1P_addition
g1_p g1_p_add(const g1_p *a, const g1_p *b) {
  if (a->infinite) return *b;
  if (b->infinite) return *a;

  if (gf_equal(a->x, b->x)) {
    GF sum_y = gf_add(a->y, b->y);
    if (gf_equal(sum_y, gf_new(0))) {
      // Pints are inverse of each other
      return g1_p_identity();
    } else {
      // Points are the same, perform doubling
      return g1_p_double(a);
    }
  }

  GF s = gf_sub(b->y, a->y);
  GF d = gf_sub(b->x, a->x);
  GF m = gf_mul(s, gf_inv(d));

  GF m_sq = gf_mul(m, m);
  GF x_r = gf_sub(gf_sub(m_sq, a->x), b->x);
  GF y_r = gf_sub(gf_mul(m, gf_sub(a->x, x_r)), a->y);

  return g1_p_new(x_r.value, y_r.value);
}

g1_p g1_p_neg(g1_p *a) {
  if (a->infinite) return *a;
  GF neg_y = gf_neg(a->y);
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

GF g1_p_generator_subgroup_size() {
  return f101(17);
}

#endif
