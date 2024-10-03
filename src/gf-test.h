#include <assert.h>
#include <stdio.h>
#include "fe.h"

void test_base() {
  u8_fe zero = f101(0);
  u8_fe one = f101(1);
  u8_fe hundred = f101(100);
  u8_fe four = f101(4);
  u8_fe twelve = f101(12);

  assert(u8_fe_equal(f101(200), u8_fe_add(hundred, hundred)));
  assert(u8_fe_equal(f101(100), u8_fe_sub(zero, one)));
  assert(u8_fe_equal(f101(0), u8_fe_div(one, zero)));

  u8_fe r = u8_fe_mul(twelve, u8_fe_div(four, twelve));
  assert(u8_fe_equal(four, r));
}

void test_vector() {
  u8_fe one = f101(1);
  u8_fe hundred = f101(100);
  assert(u8_fe_equal(one, u8_fe_pow(hundred, 0)));

  u8_fe neg_one = u8_fe_neg(one);
  assert(u8_fe_equal(hundred, neg_one));

  u8_fe two = f101(2);
  u8_fe neg_div_r = u8_fe_neg(u8_fe_div(one, two));
  assert(u8_fe_equal(f101(50), neg_div_r));

  u8_fe five = f101(5);
  assert(u8_fe_equal(f101(20), u8_fe_neg(u8_fe_div(one, five))));

  u8_fe hundred_sq = u8_fe_mul(hundred, hundred);
  u8_fe two_sq = u8_fe_mul(hundred, hundred);
  assert(u8_fe_equal(hundred_sq, two_sq));

  u8_fe hundred_cube = u8_fe_mul(hundred_sq, hundred);
  u8_fe three_pow = u8_fe_pow(hundred, 3);
  assert(u8_fe_equal(hundred_cube, three_pow));
}

int main() {
  test_base();
  test_vector();
  return 0;
}
