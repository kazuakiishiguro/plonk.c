#include <assert.h>
#include "poly.h"

#define ASSERT_POLY(p1, p2)				 \
  do {							 \
    assert(p1.len == p2.len);				 \
    for (size_t i = 0; i < p1.len; i++)			 \
	 assert(u8_fe_equal(p1.coeffs[i], p2.coeffs[i]));	\
  } while (0)

void poly_print(const poly *p) {
     bool first = true;
     for (size_t i = 0; i < p->len; i++) {
	  u8_fe coeff = p->coeffs[i];
	  if (!u8_fe_equal(coeff, u8_fe_new(0))) {
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
  u8_fe p123[] = {f101(1), f101(2), f101(3)};
  poly a = poly_new(p123, 3);
  u8_fe p14[] = {f101(1), f101(4)};
  poly b = poly_new(p14, 2);
  poly sum = poly_add(&a, &b);
  u8_fe p263[] = {f101(2), f101(6), f101(3)};
  poly expected = poly_new(p263,3);
  ASSERT_POLY(sum, expected);

  u8_fe p12345[] = {f101(1), f101(2), f101(3), f101(4), f101(5)};
  poly c = poly_new(p12345, 5);
  sum = poly_add(&a, &c);
  u8_fe p24645[] = {f101(2), f101(4), f101(6), f101(4), f101(5)};
  expected = poly_new(p24645, 5);
  ASSERT_POLY(sum, expected);

  u8_fe p12346[] = {f101(1), f101(2), f101(3), f101(4), f101(6)};
  poly d = poly_new(p12346, 5);
  sum = poly_add(&a, &d);
  u8_fe p24646[] = {f101(2), f101(4), f101(6), f101(4), f101(6)};
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
  u8_fe p123[] = {f101(1), f101(2), f101(3)};
  poly a = poly_new(p123, 3);
  poly diff = poly_sub(&a, &a);
  poly expected = poly_zero();
  ASSERT_POLY(diff, expected);

  u8_fe p12[] = {f101(1), f101(2)};
  poly b = poly_new(p12, 2);
  diff = poly_sub(&a, &b);
  u8_fe p003[] = {f101(0), f101(0), f101(3)};
  expected = poly_new(p003, 3);
  ASSERT_POLY(diff, expected);

  poly_free(&a);
  poly_free(&b);
  poly_free(&diff);
  poly_free(&expected);
}

int main() {
  test_poly_add();
  test_poly_sub();
  return 0;
}
