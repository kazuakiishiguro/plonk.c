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

int main() {
  test_base();
  test_vector();
  return 0;
}
