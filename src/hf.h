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

HF hf_new(uint8_t value) {
  int64_t tmp = value % MODULO_HF;
  if (tmp < 0) {
    tmp += MODULO_HF;
  }

  HF hf;
  hf.value = value % MODULO_HF;
  return hf;
}

HF f17(int64_t n) {
  return hf_new(n);
}

HF hf_zero() {
  return hf_new(0);
}

HF hf_one() {
  return hf_new(1);
}

bool hf_equal(HF a, HF b) {
  return a.value == b.value;
}

HF hf_add(HF a, HF b) {
  uint16_t sum = a.value + b.value;
  if (sum >= MODULO_HF) sum -= MODULO_HF;
  HF r;
  r.value = (uint8_t)sum;
  return r;
}

HF hf_sub(HF a, HF b) {
  int16_t diff = (int16_t)a.value - (int16_t)b.value;
  if (diff < 0) diff += MODULO_HF;
  HF r;
  r.value = (uint8_t)diff;
  return r;
}

HF hf_mul(HF a, HF b) {
  uint16_t product = (uint16_t)a.value * (uint16_t)b.value;
  HF r;
  r.value = product % MODULO_HF;
  return r;
}

HF hf_neg(HF a) {
  if (a.value == 0) return hf_zero();
  HF r;
  r.value = MODULO_HF - a.value;
  return r;
}

HF hf_pow(HF field, uint64_t exp) {
  HF r = hf_new(1);
  HF base = field;
  while (exp > 0) {
    if (exp % 2 == 1) {
      r = hf_mul(r, base);
    }
    exp = exp >> 1;
    base = hf_mul(base, base);
  }
  return r;
}

HF hf_inv(HF field) {
  // Using Fermat's Little Theorem: a^(p-2) mod p
  return hf_pow(field, MODULO_HF - 2);
}

HF hf_div(HF a, HF b) {
  return hf_mul(a, hf_inv(b));
}

#endif // HF_H
