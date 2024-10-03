#include <assert.h>
#include "hf.h"

void test_base() {
  HF zero = f17(0);
  HF one = f17(1);
  HF sixteen = f17(16);
  HF four = f17(4);
  HF twelve = f17(12);

  assert(hf_equal(f17(32), hf_add(sixteen, sixteen)));
  assert(hf_equal(f17(16), hf_sub(zero, one)));
  assert(hf_equal(f17(0), hf_div(one, zero)));

  HF r = hf_mul(twelve, hf_div(four, twelve));
  assert(hf_equal(four, r));
}

void test_vector() {
  HF one = f17(1);
  HF sixteen = f17(16);
  assert(hf_equal(one, hf_pow(sixteen, 0)));

  HF neg_one = hf_neg(one);
  assert(hf_equal(sixteen, neg_one));

  HF two = hf_new(2);
  HF neg_div_r = hf_neg(hf_div(one, two));
  assert(hf_equal(f17(8), neg_div_r));

  HF sixteen_sq = hf_mul(sixteen, sixteen);
  HF two_pow = hf_pow(sixteen, 2);
  assert(hf_equal(sixteen_sq, two_pow));

  HF sixteen_cube = hf_mul(sixteen_sq, sixteen);
  HF three_pow = hf_pow(sixteen, 3);
  assert(hf_equal(sixteen_cube, three_pow));

}

int main() {
  test_base();
  test_vector();
  return 0;
}
