#include <assert.h>
#include <stdio.h>
#include "fe.h"

void test_base() {
  u64_fe zero = f101(0);
  u64_fe one = f101(1);
  u64_fe hundred = f101(100);
  u64_fe four = f101(4);
  u64_fe twelve = f101(12);

  assert(f101(200).value == u64_fe_add(&hundred, &hundred).value);
  assert(f101(100).value == u64_fe_sub(&zero, &one).value);
  assert(f101(0).value == u64_fe_div(&one, &zero).value);

  u64_fe four_div_twelve = u64_fe_div(&four, &twelve);
  u64_fe r = u64_fe_mul(&twelve, &four_div_twelve);
  assert(four.value == r.value);
}

void test_vector() {
  u64_fe one = f101(1);
  u64_fe hundred = f101(100);

  u64_fe neg_one = u64_fe_sub_assign(&one);
  assert(hundred.value == neg_one.value);

  u64_fe two = f101(2);
  u64_fe div_result = u64_fe_div(&one, &two);
  u64_fe neg_div_result = u64_fe_sub_assign(&div_result);
  assert(f101(50).value == neg_div_result.value);

  u64_fe five = f101(5);
  div_result = u64_fe_div(&one, &five);
  neg_div_result = u64_fe_sub_assign(&div_result);
  assert(f101(20).value == neg_div_result.value);

  u64_fe zero_pow = u64_fe_pow(&hundred, 0);
  assert(one.value == zero_pow.value);

  u64_fe hundred_sq = u64_fe_mul(&hundred, &hundred);
  u64_fe two_pow = u64_fe_pow(&hundred, 2);
  assert(hundred_sq.value == two_pow.value);

  u64_fe hundred_cube = u64_fe_mul(&hundred_sq, &hundred);
  u64_fe three_pow = u64_fe_pow(&hundred, 3);
  assert(hundred_cube.value == three_pow.value);
}

int main() {
  test_base();
  test_vector();
  return 0;
}
