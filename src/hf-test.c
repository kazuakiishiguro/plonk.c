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

void test_hf_add() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 0; b < 17; b++) {
      HF r = hf_add(hf_new(a), hf_new(b));
      HF expected = hf_new(a+b);
      assert(hf_equal(r, expected));
    }
  }
}

void test_hf_sub() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 0; b < 17; b++) {
      HF r = hf_sub(hf_new(a), hf_new(b));
      HF expected = hf_new(a-b);
      assert(hf_equal(r, expected));
    }
  }
}

void test_hf_mul() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 0; b < 17; b++) {
      HF r = hf_mul(hf_new(a), hf_new(b));
      HF expected = hf_new(a * b);
      assert(hf_equal(r, expected));
    }
  }
}

void test_hf_neg() {
  for (uint8_t a = 0; a < 17; a++) {
    HF r = hf_neg(hf_new(a));
    uint8_t expected = (a == 0) ? 0 : (17 - a);
    assert(r.value == expected);
  }
}

void test_hf_pow() {
  for (uint8_t base = 0; base < 17; base++) {
    for (uint8_t exp = 0; exp < 17; exp++) {
      HF r = hf_pow(hf_new(base), exp);
      uint8_t expected = 1;
      for (uint8_t i = 0; i < exp; i++) {
        expected = (expected * base) % 17;
      }
      assert(r.value == expected);
    }
  }
}

void test_hf_inv() {
  for (uint8_t a = 1; a < 17; a++) { // skip zero as it doesn't have an inverse
    HF r = hf_inv(hf_new(a));
    uint8_t expected = 0;
    for (uint8_t x = 1; x < 17; x++) {
      if ((a * x) % 17 == 1) {
        expected = x;
        break;
      }
    }
    assert(r.value == expected);
  }
}

void test_hf_div() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 1; b < 17; b++) { // skip division by zero
      HF r = hf_div(hf_new(a), hf_new(b));
      uint8_t inv_b = 0;
      for (uint8_t x = 1; x < 17; x++) {
        if ((b * x) % 17 == 1) {
          inv_b = x;
          break;
        }
      }
      uint8_t expected = (a * inv_b) % 17;
      assert(r.value == expected);
    }
  }
}

int main() {
  test_base();
  test_vector();
  test_hf_add();
  test_hf_sub();
  test_hf_mul();
  test_hf_neg();
  test_hf_pow();
  test_hf_inv();
  test_hf_div();
  return 0;
}
