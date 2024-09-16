#include <stdint.h>
#include "fe.h"

typedef struct {
     u64_fe a, b;
} g2_p;

g2_p g2_p_new(uint64_t a, uint64_t b) {
     g2_p p = {
	  f101(a),
	  f101(b),
     };

     return p;
}

g2_p g2_p_generator() {
     return g2_p_new(36, 31);
}

uint64_t g2_p_embedding_degree() {
     return 2;
}

g2_p g2_p_neg(g2_p *p) {
  u64_fe neg_y = u64_fe_neg(&p ->b);
  return g2_p_new(p->a.value, neg_y.value);
}

g2_p g2_p_add(const g2_p *base, const g2_p *rhs) {
  u64_fe x, y;
  if (&base->a == &rhs->a && &base->b == &rhs->b) {
    u64_fe two = f101(2);
    u64_fe three = f101(3);
    u64_fe a_x_pow = u64_fe_pow(&base->a, 2);
    u64_fe num = u64_fe_mul(&three, &a_x_pow);
    u64_fe div = u64_fe_mul(&two, &base->b);
    u64_fe m_u = u64_fe_div(&num, &div);
    u64_fe m_u_pow = u64_fe_pow(&m_u, 2);
    u64_fe neg_two = u64_fe_neg(&two);
    u64_fe u_pow_inv = u64_fe_inv(&neg_two);
    u64_fe m_pow_2 = u64_fe_mul(&m_u_pow, &u_pow_inv);
    u64_fe two_mul_a_x = u64_fe_mul(&two, &base->a);
    x = u64_fe_sub(&m_pow_2, &two_mul_a_x);
    u64_fe inner = u64_fe_mul(&three, &base->a);
    u64_fe inner_sub_m_pow_2 = u64_fe_sub(&inner, &m_pow_2);
    u64_fe tmp = u64_fe_mul(&u_pow_inv, &m_u);
    tmp = u64_fe_mul(&tmp, &inner_sub_m_pow_2);
    y = u64_fe_sub(&tmp, &base->b);
  } else {
    u64_fe num = u64_fe_sub(&rhs->b, &base->b);
    u64_fe div = u64_fe_sub(&rhs->a, &base->a);
    u64_fe lambda_u = u64_fe_div(&num, &div);
    u64_fe two = f101(2);
    u64_fe neg_two = u64_fe_neg(&two);
    u64_fe lambda_u_pow = u64_fe_pow(&lambda_u, 2);
    lambda_u_pow = u64_fe_mul(&lambda_u_pow, &neg_two);
    x = u64_fe_sub(&lambda_u_pow, &base->a);
    x = u64_fe_sub(&x, &rhs->a);
    y = u64_fe_sub(&base->a, &x);
    y = u64_fe_mul(&lambda_u, &y);
    y = u64_fe_sub(&y, &base->b);
  }
  return g2_p_new(x.value, y.value);
}

g2_p g2_p_mul(g2_p base, u64_fe rhs) {
     uint64_t val = rhs.value;
     int flag = 0;
     g2_p result;
     while (val > 0) {
       if (val % 2 == 1) {
	 if (flag) {
	   result = g2_p_add(&result, &base);
	 } else {
	   result = base;
	   flag = 1;
	 }
       }
       val >>= 1;
       base = g2_p_add(&base, &base);
     }
     return result;
}
