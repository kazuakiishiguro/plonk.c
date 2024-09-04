#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include "utils.h"

typedef struct {
     u64_fe x, y;
     bool infinite;
} g1_p;

g1_p g1_p_new(uint64_t x_val, uint64_t y_val) {
     g1_p point = {
	  f101(x_val),
	  f101(y_val),
	  false,
     };

     return point;
}

g1_p g1_p_generator() {
     return g1_p_new(1, 2);
}

/* y^2 = x^3 + 3 */
bool g1_p_in_curve(const g1_p *point) {
     u64_fe x_cube = u64_fe_pow(&point->x, 3);
     u64_fe b = f101(3);
     return u64_fe_pow(&point->y, 2).value == u64_fe_add(&x_cube, &b).value;
}

g1_p g1_p_identity() {
     g1_p identity;
     identity.infinite = true;
     return identity;
}

g1_p g1_p_double(g1_p a) {
     u64_fe two = f101(2);
     u64_fe x_pow = u64_fe_pow(&a.x, 2);
     u64_fe three = f101(3);
     u64_fe num = u64_fe_mul(&three, &x_pow);
     u64_fe den = u64_fe_mul(&two, &a.y);
     u64_fe m = u64_fe_div(&num, &den);
     u64_fe m_pow = u64_fe_pow(&m, 2);
     u64_fe x = u64_fe_mul(&two, &a.x);
     x = u64_fe_sub(&m_pow, &x);
     u64_fe y = u64_fe_mul(&three, &a.x);
     y = u64_fe_sub(&y, &m_pow);
     y = u64_fe_mul(&m, &y);
     y = u64_fe_sub(&y, &a.y);
     g1_p r = {.x = x, .y = y, .infinite = false};
     return r;
}

g1_p g1_p_add_diff(g1_p a, g1_p b) {
     u64_fe num = u64_fe_sub(&b.y, &a.y);
     u64_fe den = u64_fe_sub(&b.x, &a.x);
     u64_fe lambda = u64_fe_div(&num, &den);
     u64_fe x = u64_fe_pow(&lambda, 2);
     x = u64_fe_sub(&x, &a.x);
     x = u64_fe_sub(&x, &b.x);
     u64_fe y = u64_fe_sub(&a.x, &x);
     y = u64_fe_mul(&lambda, &y);
     y = u64_fe_sub(&y, &a.y);
     g1_p r = {.x = x, .y = y, .infinite = false};
     return r;
}

// https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#G1P_addition
g1_p g1_p_add(g1_p a, g1_p b) {
     if (a.infinite) return b;
     if (b.infinite) return a;
     g1_p neg_b = g1_p_new(b.x.value, u64_fe_sub_assign(&b.y).value);
     if (a.x.value == neg_b.x.value && a.y.value == neg_b.y.value) return g1_p_identity();

     if (a.x.value == b.x.value && a.y.value == b.y.value) {
       return g1_p_double(a);
     } else {
       return g1_p_add_diff(a, b);
     }
}

