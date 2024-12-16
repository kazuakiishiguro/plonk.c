#ifndef HF_H
#define HF_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#define MODULO_HF 17

// Field element in HF (mod 17)
typedef struct {
  uint8_t value;
} HF;

static inline HF hf_new(int64_t value) {
  int64_t tmp = value % MODULO_HF;
  if (tmp < 0) {
    tmp += MODULO_HF;
  }

  HF hf;
  hf.value = (uint8_t)tmp;
  return hf;
}

static inline HF f17(int64_t n) {
  return hf_new(n);
}

static inline HF hf_zero() {
  HF r = {0};
  return r;
}

static inline HF hf_one() {
  HF r = {1};
  return r;
}

static inline bool hf_equal(HF a, HF b) {
  return a.value == b.value;
}

static inline HF hf_add(HF a, HF b) {
  uint8_t sum = a.value + b.value;
  if (sum >= MODULO_HF) sum -= MODULO_HF;
  HF r = { sum };
  return r;
}

static inline HF hf_sub(HF a, HF b) {
  int8_t diff = (int8_t)a.value - (int8_t)b.value;
  if (diff < 0) diff += MODULO_HF;
  HF r = { (uint8_t)diff };
  return r;
}

static inline HF hf_mul(HF a, HF b) {
  uint8_t product = (uint8_t)((uint16_t)a.value * (uint16_t)b.value % MODULO_HF);
  HF r = { product };
  return r;
}

static inline HF hf_neg(HF a) {
  HF r = { a.value == 0 ? 0 : (uint8_t)(MODULO_HF - a.value) };
  return r;
}

static inline HF hf_pow(HF base, uint64_t exp) {
  HF r = hf_one();
  while (exp > 0) {
    if (exp & 1) {
      r = hf_mul(r, base);
    }
    base = hf_mul(base, base);
    exp >>= 1;
  }
  return r;
}

// Lookup table for MODULO_HF
// This turns hf_inv() into an O(1) lookup rather than doing exponentiation.
static const uint8_t hf_inverses[17] = {
  0, 1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16
};

static inline HF hf_inv(HF field) {
  HF r = { hf_inverses[field.value] };
  return r;
}

static inline HF hf_div(HF a, HF b) {
  return hf_mul(a, hf_inv(b));
}

#endif // HF_H
