#ifndef GT_H
#define GT_H

#include <stdint.h>
#include "fe.h"

typedef struct {
  u8_fe a, b;
} gtp;

gtp gtp_new(u8_fe a, u8_fe b) {
  gtp p = {
    a,
    b
  };
  return p;
}

gtp gtp_neg(gtp *p) {
  return gtp_new(p->a, u8_fe_neg(p->b));
}

gtp gtp_mul(gtp *base, gtp *rhs) {
  u8_fe a = u8_fe_sub(u8_fe_mul(base->a, rhs->a), u8_fe_mul(u8_fe_mul(f101(2), base->b), rhs->b));
  u8_fe b = u8_fe_add(u8_fe_mul(base->a, rhs->b), u8_fe_mul(base->b, rhs->a));

  return gtp_new(a, b);
}

gtp gtp_pow(gtp *base, uint64_t exp) {
  gtp p;
  if (exp >= 101) {
    gtp tmp = gtp_pow(base, exp / 101);
    p = gtp_neg(&tmp);
    exp %= 101;
  } else {
    p = gtp_new(f101(1), f101(0));
  }

  gtp cur = *base;

  // mongomery reduction
  while (exp > 0) {
    if (exp % 2 == 1)
      p = gtp_mul(&p, &cur);
    exp >>= 1;
    cur = gtp_mul(&cur, &cur);
  }

  return p;
}

#endif
