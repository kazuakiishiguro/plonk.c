#include <stdbool.h>
#include <stdint.h>

#define MODULO 101

typedef struct {
  uint8_t value;
} u8_fe;

u8_fe u8_fe_new(uint8_t value) {
  u8_fe fe;
  fe.value = value % MODULO;
  return fe;
}

u8_fe f101(uint64_t n) {
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
  u8_fe r;
  r.value = a.value + b.value;
  if (r.value >= MODULO) r.value -= MODULO;
  return r;
}

u8_fe u8_fe_sub(u8_fe a, u8_fe b) {
  u8_fe r;
  r.value = a.value + MODULO - b.value;
  if (r.value >= MODULO) r.value -= MODULO;
  return r;
}

u8_fe u8_fe_mul(u8_fe a, u8_fe b) {
  u8_fe r;
  r.value = (a.value * b.value) % MODULO;
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
