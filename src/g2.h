#include <stdint.h>
#include "fe.h"

typedef struct {
     u64_fe x, y;
} g2_p;

g2_p g2_p_new(uint64_t x_val, uint64_t y_val) {
     g2_p point = {
	  f101(x_val),
	  f101(y_val),
     };

     return point;
}

g2_p g2_p_generator() {
     return g2_p_new(36, 31);
}

uint64_t g2_p_embedding_degree() {
     return 2;
}

g2_p g2_p_sub_assign(g2_p *a) {
  u64_fe neg_y = u64_fe_sub_assign(&a ->y);
  return g2_p_new(a->x.value, neg_y.value);
}

g2_p g2_p_add(const g2_p *a, const g2_p *b) {
  u64_fe x, y;
  if (&a->x == &b->x && &a->y == &b->y) {
    u64_fe two = f101(2);
    u64_fe three = f101(3);
    u64_fe a_x_pow = u64_fe_pow(&a->x, 2);
    u64_fe num = u64_fe_mul(&three, &a_x_pow);
    u64_fe div = u64_fe_mul(&two, &a->y);
    u64_fe m_u = u64_fe_div(&num, &div);
    u64_fe m_u_pow = u64_fe_pow(&m_u, 2);
    u64_fe neg_two = u64_fe_sub_assign(&two);
    u64_fe u_pow_inv = u64_fe_inv(&neg_two);
    u64_fe m_pow_2 = u64_fe_mul(&m_u_pow, &u_pow_inv);
    u64_fe two_mul_a_x = u64_fe_mul(&two, &a->x);
    x = u64_fe_sub(&m_pow_2, &two_mul_a_x);
    u64_fe inner = u64_fe_mul(&three, &a->x);
    u64_fe inner_sub_m_pow_2 = u64_fe_sub(&inner, &m_pow_2);
    u64_fe tmp = u64_fe_mul(&u_pow_inv, &m_u);
    tmp = u64_fe_mul(&tmp, &inner_sub_m_pow_2);
    y = u64_fe_sub(&tmp, &a->y);
  } else {
    u64_fe num = u64_fe_sub(&b->y, &a->y);
    u64_fe div = u64_fe_sub(&b->x, &a->x);
    u64_fe lambda_u = u64_fe_div(&num, &div);
    u64_fe two = f101(2);
    u64_fe neg_two = u64_fe_sub_assign(&two);
    u64_fe lambda_u_pow = u64_fe_pow(&lambda_u, 2);
    lambda_u_pow = u64_fe_mul(&lambda_u_pow, &neg_two);
    x = u64_fe_sub(&lambda_u_pow, &a->x);
    x = u64_fe_sub(&x, &b->x);
    y = u64_fe_sub(&a->x, &x);
    y = u64_fe_mul(&lambda_u, &y);
    y = u64_fe_sub(&y, &a->y);
  }
  return g2_p_new(x.value, y.value);
}
