#include <stdint.h>
#include "fe.h"

typedef struct {
     u64_fe a, b;
} gtp;

gtp gtp_new(u64_fe a, u64_fe b) {
     gtp p = {
	  a,
	  b
     };
     return p;
}

gtp gtp_sub_assign(gtp *p) {
     u64_fe neg_y = u64_fe_sub_assign(&p->b);
     return gtp_new(p->a, neg_y);
}

gtp gtp_mul(gtp *base, gtp *rhs) {
    u64_fe base_mul_rhs = u64_fe_mul(&base->a, &rhs->a);
    u64_fe two = f101(2);
    u64_fe two_mul_base = u64_fe_mul(&two, &base->b);
    u64_fe two_mul_base_rhs = u64_fe_mul(&two_mul_base, &rhs->b);
    u64_fe a = u64_fe_sub(&base_mul_rhs, &two_mul_base_rhs);

    u64_fe base_mul_rhs_b = u64_fe_mul(&base->a, &rhs->b);
    u64_fe base_mul_rhs_a = u64_fe_mul(&base->b, &rhs->a);
    u64_fe b = u64_fe_add(&base_mul_rhs_b, &base_mul_rhs_a);
    return gtp_new(a, b);
}
