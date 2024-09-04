#include <stdbool.h>
#include <stdint.h>

#define MODULO 101

typedef struct {
     uint64_t value;
     uint64_t modulo;
} u64_fe;

u64_fe u64_fe_new(uint64_t value, uint64_t modulo) {
     u64_fe fe;
     fe.value = value % modulo;
     fe.modulo = modulo;
     return fe;
}

u64_fe f101(uint64_t n) {
     return u64_fe_new(n, MODULO);
}

u64_fe u64_fe_one(uint64_t modulo) {
     return u64_fe_new(1, modulo);
}

static int64_t xgcd(int64_t a, int64_t b, int64_t *x, int64_t *y) {
     if (a == 0) {
	  *x = 0;
	  *y = 1;
	  return b;
     }

     int64_t x1, y1;
     int64_t gcd = xgcd(b % a, a, &x1, &y1);

     *x = y1 - (b/a) * x1;
     *y = x1;

     return gcd;
}

u64_fe u64_fe_inv(const u64_fe *field) {
     int64_t x, y;
     int64_t g = xgcd((int64_t)field->value, (int64_t)field->modulo, &x, &y);
     if (g != 1) {
	  return u64_fe_new(0, field->modulo); // no inverse operation
     } else {
	  int64_t r = (x % (int64_t)field->modulo + (int64_t)field->modulo) % (int64_t)field->modulo;
	  return u64_fe_new((uint64_t)r, field->modulo);
     }
}

inline static bool is_odd(uint64_t n) {
     return n & 1;
}

u64_fe u64_fe_pow(const u64_fe *field, uint64_t exp) {
     u64_fe r = u64_fe_one(field->modulo);
     u64_fe base = *field;
     while (exp > 0) {
	  if (is_odd(exp)) {
	       r.value = (r.value * base.value) % r.modulo;
	  }
	  exp = exp >> 1;
	  base.value = (base.value * base.value) % base.modulo;
     }
     return r;
}

u64_fe u64_fe_add(const u64_fe *a, const u64_fe *b) {
     return u64_fe_new((a->value + b->value) % a->modulo, a->modulo);
}

u64_fe u64_fe_sub(const u64_fe *a, const u64_fe *b) {
     return u64_fe_new((a->value + a->modulo - b->value) % a->modulo, a->modulo);
}

u64_fe u64_fe_sub_assign(u64_fe *a) {
     return u64_fe_new(a->modulo - a->value, a->modulo);
}

u64_fe u64_fe_mul(const u64_fe *a, const u64_fe *b) {
     return u64_fe_new((a->value * b->value) % a->modulo, a->modulo);
}

u64_fe u64_fe_div(const u64_fe *a, const u64_fe *b) {
     u64_fe inv_b = u64_fe_inv(b);
     if (inv_b.value == 0) {
	  return u64_fe_new(0, a->modulo); // no division possible
     }
     return u64_fe_mul(a, &inv_b);
}
