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

HF hf_new(uint8_t x) {
  HF result;
  result.value = x % MODULO_HF;
  return result;
}

HF f17(int64_t n) {
  return hf_new(n);
}

bool hf_equal(HF a, HF b) {
  return a.value == b.value;
}

HF hf_add(HF a, HF b) {
  return hf_new((a.value + b.value) % MODULO_HF);
}

HF hf_sub(HF a, HF b) {
  return hf_new((a.value + MODULO_HF - b.value) % MODULO_HF);
}

HF hf_mul(HF a, HF b) {
  return hf_new((a.value * b.value) % MODULO_HF);
}

HF hf_pow(HF base, uint32_t exponent) {
  HF result = hf_new(1);
  HF b = base;
  while (exponent > 0) {
    if (exponent % 2 == 1) {
      result = hf_mul(result, b);
    }
    b = hf_mul(b, b);
    exponent /= 2;
  }
  return result;
}

HF hf_inv(HF a) {
  if (a.value == 0) {
    fprintf(stderr, "Attempt to invert zero in finite field\n");
    exit(EXIT_FAILURE);
  }
  // Using Fermat's Little Theorem: a^(p-2) mod p
  return hf_pow(a, MODULO_HF - 2);
}

HF hf_div(HF a, HF b) {
  return hf_mul(a, hf_inv(b));
}

#endif // HF_H
