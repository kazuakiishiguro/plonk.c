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

int main() {
  test_base();
  return 0;
}
