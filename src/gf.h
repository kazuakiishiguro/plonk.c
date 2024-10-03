#ifndef FE_H
#define FE_H

#include <stdbool.h>
#include <stdint.h>
#include "hf.h"

#define MODULO 101

typedef struct {
  uint8_t value;
} GF;

GF gf_new(int64_t value) {
  int64_t tmp = value % MODULO;
  if (tmp < 0) {
    tmp += MODULO;
  }

  GF gf;
  gf.value = (uint8_t)tmp;
  return gf;
}

GF f101(int64_t n) {
  return gf_new(n);
}

GF gf_zero() {
  return gf_new(0);
}

GF gf_one() {
  return gf_new(1);
}

bool gf_in_field(GF a) {
  return a.value < MODULO;
}

inline static bool is_odd(uint64_t n) {
  return n & 1;
}

bool gf_equal(GF a, GF b) {
  return a.value == b.value;
}

GF gf_add(GF a, GF b) {
  uint16_t sum = a.value + b.value;
  if (sum >= MODULO) sum -= MODULO;
  GF r;
  r.value = (uint8_t)sum;
  return r;
}

GF gf_sub(GF a, GF b) {
  int16_t diff = (int16_t)a.value - (int16_t)b.value;
  if (diff < 0) diff += MODULO;
  GF r;
  r.value = (uint8_t)diff;
  return r;
}

GF gf_mul(GF a, GF b) {
  uint16_t product = (uint16_t)a.value * (uint16_t)b.value;
  GF r;
  r.value = product % MODULO;
  return r;
}

GF gf_neg(GF a) {
  if (a.value == 0) return gf_new(0);
  GF r;
  r.value = MODULO - a.value;
  return r;
}

GF gf_pow(GF field, uint64_t exp) {
  GF r = gf_one();
  GF base = field;
  while (exp > 0) {
    if (is_odd(exp)) {
      r = gf_mul(r, base);
    }
    exp = exp >> 1;
    base = gf_mul(base, base);
  }
  return r;
}

GF gf_inv(GF field) {
  // fermat's little theorem for inverse calculation
  return gf_pow(field, MODULO -2);
}

GF gf_div(GF a, GF b) {
  GF inv_b = gf_inv(b);
  return gf_mul(a, inv_b);
}

GF gf_from_hf(HF hf_elem) {
  return gf_new(hf_elem.value);
}

#endif
