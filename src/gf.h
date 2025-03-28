#ifndef FE_H
#define FE_H

#include <stdbool.h>
#include <stdint.h>
#include "hf.h"

#define MODULO_GF 101

typedef struct {
  uint8_t value;
} GF;

static inline GF gf_new(int64_t value) {
  int64_t tmp = value % MODULO_GF;
  if (tmp < 0) {
    tmp += MODULO_GF;
  }

  GF gf;
  gf.value = (uint8_t)tmp;
  return gf;
}

static inline GF f101(int64_t n) {
  return gf_new(n);
}

static inline GF gf_zero() {
  GF r = {0};
  return r;
}

static inline bool is_odd(uint64_t n) {
  return n & 1;
}

static inline GF gf_one() {
  GF r = {1};
  return r;
}

static inline bool gf_equal(GF a, GF b) {
  return a.value == b.value;
}

static inline GF gf_add(GF a, GF b) {
  uint16_t sum = a.value + b.value;
  if (sum >= MODULO_GF) sum -= MODULO_GF;
  GF r;
  r.value = (uint8_t)sum; // This cast is safe because sum < 2*MODULO_GF
  return r;
}

static inline GF gf_sub(GF a, GF b) {
  int16_t diff = (int16_t)a.value - (int16_t)b.value;
  if (diff < 0) diff += MODULO_GF;
  GF r;
  r.value = (uint8_t)diff; // Safe cast as 0 <= diff < MODULO_GF
  return r;
}

static inline GF gf_mul(GF a, GF b) {
  uint16_t product = (uint16_t)a.value * (uint16_t)b.value;
  GF r;
  r.value = product % MODULO_GF;
  return r;
}

static inline GF gf_neg(GF a) {
  if (a.value == 0) return gf_zero();
  GF r;
  r.value = a.value == 0 ? 0 : MODULO_GF - a.value;
  return r;
}

static inline GF gf_pow(GF field, uint64_t exp) {
  GF r = gf_one();
  GF base = field;
  while (exp > 0) {
    if (is_odd(exp)) {
      r = gf_mul(r, base);
    }
    exp >>= 1;
    base = gf_mul(base, base);
  }
  return r;
}

static inline GF gf_inv(GF field) {
  // Using Fermat's Little Theorem: a^(p-2) mod p
  return gf_pow(field, MODULO_GF -2);
}

static inline GF gf_div(GF a, GF b) {
  return gf_mul(a, gf_inv(b));
}

static inline GF gf_from_hf(HF hf_elem) {
  return gf_new(hf_elem.value);
}

#endif
