#include <stdbool.h>
#include <stdint.h>

#define MODULO 101

typedef struct {
  uint8_t value;
} u64_fe;

u64_fe u64_fe_new(uint64_t value) {
  u64_fe fe;
  fe.value = value % MODULO;
  return fe;
}

u64_fe f101(uint64_t n) {
  return u64_fe_new(n);
}

u64_fe u64_fe_one() {
  return u64_fe_new(1);
}

inline static bool is_odd(uint64_t n) {
  return n & 1;
}

u64_fe u64_fe_pow(const u64_fe *field, uint64_t exp) {
  u64_fe r = u64_fe_one();
  u64_fe base = *field;
  while (exp > 0) {
    if (is_odd(exp)) {
      r.value = (r.value * base.value) % MODULO;
    }
    exp = exp >> 1;
    base.value = (base.value * base.value) % MODULO;
  }
  return r;
}

u64_fe u64_fe_inv(const u64_fe *field) {
  // fermat's little theorem for inverse calculation
  return u64_fe_pow(field, MODULO -2);
}

u64_fe u64_fe_add(const u64_fe *a, const u64_fe *b) {
  u64_fe r;
  r.value = a->value + b->value;
  if (r.value >= MODULO) r.value -= MODULO;
  return r;
}

u64_fe u64_fe_sub(const u64_fe *a, const u64_fe *b) {
  u64_fe r;
  r.value = a->value + MODULO - b->value;
  if (r.value >= MODULO) r.value -= MODULO;
  return r;
}

u64_fe u64_fe_mul(const u64_fe *a, const u64_fe *b) {
  u64_fe r;
  r.value = (a->value * b->value) % MODULO;
  return r;
}

u64_fe u64_fe_neg(const u64_fe *a) {
  u64_fe r;
  r.value = (MODULO - a->value) % MODULO;
  return r;
}

u64_fe u64_fe_div(const u64_fe *a, const u64_fe *b) {
  u64_fe inv_b = u64_fe_inv(b);
  return u64_fe_mul(a, &inv_b);
}
