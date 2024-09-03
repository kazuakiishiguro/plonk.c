#include <assert.h>
#include <stdio.h>
#include "utils.h"

#define MODULO 101

u64_fe f101(uint64_t n) {
  return u64_fe_new(n, MODULO);
}

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

int main() {
  test_base();
  return 0;
}
