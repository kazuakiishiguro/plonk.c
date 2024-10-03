#ifndef GT_H
#define GT_H

#include <stdint.h>
#include "gf.h"

typedef struct {
  GF a, b;
} GTP;

GTP gtp_new(GF a, GF b) {
  GTP p = {
    a,
    b
  };
  return p;
}

GTP gtp_neg(GTP *p) {
  return gtp_new(p->a, gf_neg(p->b));
}

GTP gtp_mul(GTP *base, GTP *rhs) {
  GF a = gf_sub(gf_mul(base->a, rhs->a), gf_mul(gf_mul(f101(2), base->b), rhs->b));
  GF b = gf_add(gf_mul(base->a, rhs->b), gf_mul(base->b, rhs->a));

  return gtp_new(a, b);
}

GTP gtp_pow(GTP *base, uint64_t exp) {
  GTP p;
  if (exp >= 101) {
    GTP tmp = gtp_pow(base, exp / 101);
    p = gtp_neg(&tmp);
    exp %= 101;
  } else {
    p = gtp_new(f101(1), f101(0));
  }

  GTP cur = *base;

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
