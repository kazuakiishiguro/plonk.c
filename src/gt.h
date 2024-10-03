#ifndef GT_H
#define GT_H

#include <stdint.h>
#include "gf.h"

typedef struct {
  GF a, b;
} gtp;

gtp gtp_new(GF a, GF b) {
  gtp p = {
    a,
    b
  };
  return p;
}

gtp gtp_neg(gtp *p) {
  return gtp_new(p->a, gf_neg(p->b));
}

gtp gtp_mul(gtp *base, gtp *rhs) {
  GF a = gf_sub(gf_mul(base->a, rhs->a), gf_mul(gf_mul(f101(2), base->b), rhs->b));
  GF b = gf_add(gf_mul(base->a, rhs->b), gf_mul(base->b, rhs->a));

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
