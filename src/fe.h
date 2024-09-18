#ifndef FE_H
#define FE_H

#include <stdbool.h>
#include <stdint.h>

#define MODULO 101

typedef struct {
  uint8_t value;
} u8_fe;

u8_fe u8_fe_new(int64_t value) {
  int64_t tmp = value % MODULO;
  if (tmp < 0) {
    tmp += MODULO;
  }

  u8_fe fe;
  fe.value = (uint8_t)tmp;
  return fe;
}

u8_fe f101(int64_t n) {
  return u8_fe_new(n);
}

u8_fe u8_fe_one() {
  return u8_fe_new(1);
}

inline static bool is_odd(uint64_t n) {
  return n & 1;
}

bool u8_fe_equal(u8_fe a, u8_fe b) {
  return a.value == b.value;
}

u8_fe u8_fe_add(u8_fe a, u8_fe b) {
  uint16_t sum = a.value + b.value;
  if (sum >= MODULO) sum -= MODULO;
  u8_fe r;
  r.value = (uint8_t)sum;
  return r;
}

u8_fe u8_fe_sub(u8_fe a, u8_fe b) {
  int16_t diff = (int16_t)a.value - (int16_t)b.value;
  if (diff < 0) diff += MODULO;
  u8_fe r;
  r.value = (uint8_t)diff;
  return r;
}

u8_fe u8_fe_mul(u8_fe a, u8_fe b) {
  uint16_t product = (uint16_t)a.value * (uint16_t)b.value;
  u8_fe r;
  r.value = product % MODULO;
  return r;
}

u8_fe u8_fe_neg(u8_fe a) {
  u8_fe r;
  r.value = (MODULO - a.value) % MODULO;
  return r;
}

u8_fe u8_fe_pow(u8_fe field, uint64_t exp) {
  u8_fe r = u8_fe_one();
  u8_fe base = field;
  while (exp > 0) {
    if (is_odd(exp)) {
      r.value = (r.value * base.value) % MODULO;
    }
    exp = exp >> 1;
    base.value = (base.value * base.value) % MODULO;
  }
  return r;
}

u8_fe u8_fe_inv(u8_fe field) {
  // fermat's little theorem for inverse calculation
  return u8_fe_pow(field, MODULO -2);
}

u8_fe u8_fe_div(u8_fe a, u8_fe b) {
  u8_fe inv_b = u8_fe_inv(b);
  return u8_fe_mul(a, inv_b);
}

#endif
