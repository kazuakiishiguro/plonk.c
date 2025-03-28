#include <assert.h>
#include <stdio.h>
#include "gf.h"

void test_base() {
  GF zero = f101(0);
  GF one = f101(1);
  GF hundred = f101(100);
  GF four = f101(4);
  GF twelve = f101(12);

  assert(gf_equal(f101(200), gf_add(hundred, hundred)));
  assert(gf_equal(f101(100), gf_sub(zero, one)));
  assert(gf_equal(f101(0), gf_div(one, zero)));

  GF r = gf_mul(twelve, gf_div(four, twelve));
  assert(gf_equal(four, r));
}

void test_vector() {
  GF one = f101(1);
  GF hundred = f101(100);
  assert(gf_equal(one, gf_pow(hundred, 0)));

  GF neg_one = gf_neg(one);
  assert(gf_equal(hundred, neg_one));

  GF two = f101(2);
  GF neg_div_r = gf_neg(gf_div(one, two));
  assert(gf_equal(f101(50), neg_div_r));

  GF five = f101(5);
  assert(gf_equal(f101(20), gf_neg(gf_div(one, five))));

  GF hundred_sq = gf_mul(hundred, hundred);
  GF two_pow = gf_pow(hundred, 2);
  assert(gf_equal(hundred_sq, two_pow));

  GF hundred_cube = gf_mul(hundred_sq, hundred);
  GF three_pow = gf_pow(hundred, 3);
  assert(gf_equal(hundred_cube, three_pow));
}

void test_corners() {
  GF zero = gf_zero();
  GF one = gf_one();
  GF two = f101(2);
  GF five = f101(5);
  GF last = f101(MODULO_GF - 1); // 100

  // Zero interactions
  assert(gf_equal(five, gf_add(five, zero)));
  assert(gf_equal(five, gf_sub(five, zero)));
  assert(gf_equal(zero, gf_sub(five, five)));
  assert(gf_equal(zero, gf_mul(five, zero)));
  assert(gf_equal(zero, gf_mul(zero, five)));
  assert(gf_equal(zero, gf_div(zero, five))); // 0 / non-zero
  assert(gf_equal(zero, gf_neg(zero)));
  assert(gf_equal(zero, gf_pow(zero, 5))); // 0^k for k>0
  assert(gf_equal(one, gf_pow(zero, 0))); // 0^0 convention/implementation = 1

  // One interactions
  assert(gf_equal(five, gf_mul(five, one)));
  assert(gf_equal(five, gf_mul(one, five)));
  assert(gf_equal(five, gf_div(five, one)));
  assert(gf_equal(one, gf_div(five, five)));
  assert(gf_equal(one, gf_pow(one, 5)));
  assert(gf_equal(five, gf_pow(five, 1)));
  assert(gf_equal(one, gf_inv(one)));

  // Modulus boundaries
  assert(gf_equal(zero, gf_add(last, one)));
  assert(gf_equal(last, gf_sub(zero, one))); // Already in test_base
  assert(gf_equal(last, gf_sub(one, two)));

  // Inverse/Negation
  assert(gf_equal(one, gf_mul(five, gf_inv(five))));
  assert(gf_equal(five, gf_neg(gf_neg(five))));
  assert(gf_equal(zero, gf_add(five, gf_neg(five))));

  // gf_new
  assert(gf_equal(zero, f101(MODULO_GF)));
  assert(gf_equal(one, f101(MODULO_GF + 1)));
  assert(gf_equal(last, f101(-1)));
  assert(gf_equal(zero, f101(-MODULO_GF)));
  assert(gf_equal(f101(5), f101(1010 + 5)));
  assert(gf_equal(f101(96), f101(-1010 - 5)));
}

int main() {
  test_base();
  test_vector();
  test_corners();
  return 0;
}
