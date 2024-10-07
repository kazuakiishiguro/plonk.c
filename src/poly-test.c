#include <assert.h>
#include "gf.h"
#include "poly.h"

#define ASSERT_POLY(p1, p2)				 \
  do {							 \
    assert(p1.len == p2.len);				 \
    for (size_t i = 0; i < p1.len; i++)			 \
	 assert(hf_equal(p1.coeffs[i], p2.coeffs[i]));	\
  } while (0)

void poly_print(const POLY *p) {
     bool first = true;
     for (size_t i = 0; i < p->len; i++) {
	  HF coeff = p->coeffs[i];
	  if (!hf_equal(coeff, hf_new(0))) {
	       if (!first)
		    printf(" + ");
	       if (i == 0) {
		    printf("%u", coeff.value);
	       } else if (i == 1){
		    printf("%u*x", coeff.value);
	       } else {
		    if (coeff.value != 1)
			 printf("%u*", coeff.value);
		    printf("x^%zu", i);
	       }
	       first = false;
	  }
     }
     if (first) {
	  printf("0");
     }
     printf("\n");
}

void test_poly_add() {
  HF p123[] = {f17(1), f17(2), f17(3)};
  POLY a = poly_new(p123, 3);
  HF p14[] = {f17(1), f17(4)};
  POLY b = poly_new(p14, 2);
  POLY sum = poly_add(&a, &b);
  HF p263[] = {f17(2), f17(6), f17(3)};
  POLY expected = poly_new(p263,3);
  ASSERT_POLY(sum, expected);

  HF p12345[] = {f17(1), f17(2), f17(3), f17(4), f17(5)};
  POLY c = poly_new(p12345, 5);
  sum = poly_add(&a, &c);
  HF p24645[] = {f17(2), f17(4), f17(6), f17(4), f17(5)};
  expected = poly_new(p24645, 5);
  ASSERT_POLY(sum, expected);

  HF p12346[] = {f17(1), f17(2), f17(3), f17(4), f17(6)};
  POLY d = poly_new(p12346, 5);
  sum = poly_add(&a, &d);
  HF p24646[] = {f17(2), f17(4), f17(6), f17(4), f17(6)};
  expected = poly_new(p24646, 5);
  ASSERT_POLY(sum, expected);

  poly_free(&a);
  poly_free(&b);
  poly_free(&c);
  poly_free(&d);
  poly_free(&sum);
  poly_free(&expected);
}

void test_poly_sub() {
  HF p123[] = {f17(1), f17(2), f17(3)};
  POLY a = poly_new(p123, 3);
  POLY diff = poly_sub(&a, &a);
  POLY expected = poly_zero();
  ASSERT_POLY(diff, expected);

  HF p12[] = {f17(1), f17(2)};
  POLY b = poly_new(p12, 2);
  diff = poly_sub(&a, &b);
  HF p003[] = {f17(0), f17(0), f17(3)};
  expected = poly_new(p003, 3);
  ASSERT_POLY(diff, expected);

  poly_free(&a);
  poly_free(&b);
  poly_free(&diff);
  poly_free(&expected);
}

void test_poly_mul() {
  HF ca[] = {f17(5), f17(0), f17(10), f17(6)};
  POLY a = poly_new(ca, 4);
  HF cb[] = {f17(1), f17(2), f17(4)};
  POLY b = poly_new(cb, 3);
  POLY product = poly_mul(&a, &b);
  HF cexp[] = {f17(5), f17(10), f17(30), f17(26), f17(52), f17(24)};
  POLY expected = poly_new(cexp, 6);
  ASSERT_POLY(product, expected);

  poly_free(&a);
  poly_free(&b);
  poly_free(&product);
  poly_free(&expected);
}

void test_poly_eval() {
  HF coeffs[] = {f17(1), f17(2), f17(1)};
  POLY p = poly_new(coeffs, 3);
  HF x = hf_new(2);
  HF eval = poly_eval(&p, x);
  HF expected = f17(9);
  hf_equal(eval, expected);
}

void test_poly_z() {
  HF points[] = {f17(1), f17(5)};
  POLY z = poly_z(points, 2);
  HF c[] = {f17(5), f17(-6), f17(1)};
  POLY expected = poly_new(c, 3);
  ASSERT_POLY(z, expected);

  poly_free(&z);
  poly_free(&expected);
}

void test_poly_lagrange() {
  HF x_points[] = {f17(1), f17(5), f17(7), f17(3)};
  HF y_points[] = {f17(2), f17(7), f17(9), f17(1)};
  size_t num_points = sizeof(x_points) / sizeof(x_points[0]);
  POLY l = poly_lagrange(x_points, y_points, num_points);
  for (size_t i = 0; i < num_points;  i++) {
    HF x = x_points[i];
    HF expected_y = y_points[i];
    HF actual_y =  poly_eval(&l, x);
    assert(hf_equal(actual_y, expected_y));
  }
  poly_free(&l);
}

int main() {
  test_poly_add();
  test_poly_sub();
  test_poly_eval();
  test_poly_z();
  test_poly_lagrange();
  return 0;
}
