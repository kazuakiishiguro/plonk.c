#ifndef HF_H
#define HF_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define MODULO_HF 17

// Field element in HF (mod 17)
typedef struct {
  uint8_t value;
} hf_fe;

hf_fe hf_fe_new(uint8_t x) {
  hf_fe result;
  result.value = x % MODULO_HF;
  return result;
}

hf_fe f17(int64_t n) {
  return hf_fe_new(n);
}

hf_fe hf_fe_add(hf_fe a, hf_fe b) {
  return hf_fe_new((a.value + b.value) % MODULO_HF);
}

hf_fe hf_fe_sub(hf_fe a, hf_fe b) {
  return hf_fe_new((a.value + MODULO_HF - b.value) % MODULO_HF);
}

hf_fe hf_fe_mul(hf_fe a, hf_fe b) {
  return hf_fe_new((a.value * b.value) % MODULO_HF);
}

hf_fe hf_fe_pow(hf_fe base, uint32_t exponent) {
  hf_fe result = hf_fe_new(1);
  hf_fe b = base;
  while (exponent > 0) {
    if (exponent % 2 == 1) {
      result = hf_fe_mul(result, b);
    }
    b = hf_fe_mul(b, b);
    exponent /= 2;
  }
  return result;
}

hf_fe hf_fe_inv(hf_fe a) {
  if (a.value == 0) {
    fprintf(stderr, "Attempt to invert zero in finite field\n");
    exit(EXIT_FAILURE);
  }
  // Using Fermat's Little Theorem: a^(p-2) mod p
  return hf_fe_pow(a, MODULO_HF - 2);
}

hf_fe hf_fe_div(hf_fe a, hf_fe b) {
  return hf_fe_mul(a, hf_fe_inv(b));
}

#endif // HF_H
