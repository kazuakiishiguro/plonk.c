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

  // Division by zero is not mathematically defined in a field, but the code returns zero.
  // This test checks current behavior (not truly correct in terms of fields):
  assert(hf_equal(f17(0), hf_div(one, zero)));

  // Test a*inv(a) * b == b
  HF r = hf_mul(twelve, hf_div(four, twelve));
  assert(hf_equal(four, r));
}

void test_vector() {
  HF one = f17(1);
  HF sixteen = f17(16);
  assert(hf_equal(one, hf_pow(sixteen, 0))); // x^0 = 1  

  HF neg_one = hf_neg(one);
  assert(hf_equal(sixteen, neg_one)); // -1 = 16 mod 17

  HF two = hf_new(2);
  HF one_div_two = hf_div(one, two); // 1/2 = 9 mod 17
  assert(hf_equal(f17(9), one_div_two));
  HF neg_div_r = hf_neg(one_div_two); // -(1/2) = -9 = 8 mod 17
  assert(hf_equal(f17(8), neg_div_r));

  HF sixteen_sq = hf_mul(sixteen, sixteen); // (-1)*(-1) = 1 mod 17
  HF two_pow = hf_pow(sixteen, 2);
  assert(hf_equal(sixteen_sq, two_pow));
  assert(hf_equal(f17(1), sixteen_sq)); // Verify (-1)^2 = 1

  HF sixteen_cube = hf_mul(sixteen_sq, sixteen); // 1 * (-1) = -1 = 16 mod 17
  HF three_pow = hf_pow(sixteen, 3);
  assert(hf_equal(sixteen_cube, three_pow));
  assert(hf_equal(f17(16), sixteen_cube)); // Verify (-1)^3 = -1
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
      HF expected = hf_new(a-b); // hf_new handles the negative result correctly
      assert(hf_equal(r, expected));
    }
  }
}

void test_hf_mul() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 0; b < 17; b++) {
      HF r = hf_mul(hf_new(a), hf_new(b));
      // Use wider type for intermediate multiplication
      HF expected = hf_new((int64_t)a * b);
      assert(hf_equal(r, expected));
    }
  }
}

void test_hf_neg() {
  for (uint8_t a = 0; a < 17; a++) {
    HF r = hf_neg(hf_new(a));
    uint8_t expected_val = (a == 0) ? 0 : (17 - a);
    HF expected = hf_new(expected_val);
    assert(hf_equal(r, expected));
    // Also check definition: a + (-a) == 0
    assert(hf_equal(hf_add(hf_new(a), r), hf_zero()));
  }
}

void test_hf_pow() {
  for (uint8_t base = 0; base < 17; base++) {
    // Test small exponents exhaustively, plus some larger ones
    uint64_t exponents[] = {0, 1, 2, 5, 15, 16, 17, 30, 65};
    int num_exponents = sizeof(exponents) / sizeof(exponents[0]);

    for (int i = 0; i < num_exponents; i++) {
      uint64_t exp = exponents[i];
      HF r = hf_pow(hf_new(base), exp);

      // Calculate expected value using simpler multiplication loop
      // Use wider intermediate type for multiplication before modulo
      uint64_t expected_val_long = 1;
      for (uint64_t j = 0; j < exp; j++) {
        expected_val_long = (expected_val_long * base) % 17;
      }
      // Handle base=0, exp=0 case (0^0 = 1)
      if (base == 0 && exp == 0) {
          expected_val_long = 1;
      } else if (base == 0 && exp > 0) {
          expected_val_long = 0;
      }

      HF expected = hf_new((int64_t)expected_val_long);
      assert(hf_equal(r, expected));
    }
  }
  // Explicit check for 0^0 = 1
  assert(hf_equal(hf_pow(hf_zero(), 0), hf_one()));
  // Explicit check for 0^n = 0 for n>0
  assert(hf_equal(hf_pow(hf_zero(), 1), hf_zero()));
  assert(hf_equal(hf_pow(hf_zero(), 5), hf_zero()));
}

void test_hf_inv() {
  // Test against manual calculation of inverse
  for (uint8_t a = 1; a < 17; a++) { // skip zero as it doesn't have an inverse
    HF r = hf_inv(hf_new(a));
    uint8_t expected_val = 0;

    // Find inverse by checking all possibilities
    for (uint8_t x = 1; x < 17; x++) {
      if (((uint16_t)a * x) % 17 == 1) {
        expected_val = x;
        break;
      }
    }
    HF expected = hf_new(expected_val);
    assert(hf_equal(r, expected));
  }

  // Additional verification: (a * inv(a)) == 1 for all a != 0
  for (uint8_t a = 1; a < 17; a++) {
    HF A = hf_new(a);
    HF inv_a = hf_inv(A);
    HF check = hf_mul(A, inv_a);
    assert(hf_equal(check, hf_one()));
  }

  // Check inv(0) returns 0 (as per lookup table)
  assert(hf_equal(hf_inv(hf_zero()), hf_zero()));
}

void test_hf_div() {
  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 1; b < 17; b++) { // skip division by zero
      HF r = hf_div(hf_new(a), hf_new(b));

      // Calculate expected value a * inv(b)
      uint8_t inv_b_val = 0;
      for (uint8_t x = 1; x < 17; x++) {
        if (((uint16_t)b * x) % 17 == 1) {
          inv_b_val = x;
          break;
        }
      }
      // Use wider intermediate type for multiplication before modulo
      uint8_t expected_val = ((uint16_t)a * inv_b_val) % 17;
      HF expected = hf_new(expected_val);
      assert(hf_equal(r, expected));
    }
  }
  // Check division by zero results in zero (current behavior)
   assert(hf_equal(hf_div(hf_new(5), hf_zero()), hf_zero()));
   assert(hf_equal(hf_div(hf_one(), hf_zero()), hf_zero()));
}

void test_field_axioms() {
  HF Z = hf_zero();
  HF O = hf_one();

  for (uint8_t a = 0; a < 17; a++) {
    for (uint8_t b = 0; b < 17; b++) {
      for (uint8_t c = 0; c < 17; c++) {
        HF A = hf_new(a), B = hf_new(b), C = hf_new(c);
        // Associativity of addition: (A + B) + C == A + (B + C)
        assert(hf_equal(hf_add(hf_add(A, B), C), hf_add(A, hf_add(B, C))));
        // Associativity of multiplication: (A * B) * C == A * (B * C)
        assert(hf_equal(hf_mul(hf_mul(A, B), C), hf_mul(A, hf_mul(B, C))));
        // Commutativity of addition: A + B == B + A
        assert(hf_equal(hf_add(A, B), hf_add(B, A)));
        // Commutativity of multiplication: A * B == B * A
        assert(hf_equal(hf_mul(A, B), hf_mul(B, A)));
        // Distributivity: A * (B + C) == (A * B) + (A * C)
        assert(hf_equal(hf_mul(A, hf_add(B, C)), hf_add(hf_mul(A, B), hf_mul(A, C))));
      }
    }
  }
}

void test_identities_inverses() {
    for (uint8_t a = 0; a < 17; a++) {
        HF A = hf_new(a);
        HF Z = hf_zero();
        HF O = hf_one();

        // Additive identity: A+0 = A
        assert(hf_equal(hf_add(A, Z), A));
        // Multiplicative identity: A*1 = A
        assert(hf_equal(hf_mul(A, O), A));
        // Additive inverse property: A+(-A) = 0
        assert(hf_equal(hf_add(A, hf_neg(A)), Z));
        // Multiplication by zero: A*0 = 0
        assert(hf_equal(hf_mul(A, Z), Z));
        // Division by one: A/1 = A
        assert(hf_equal(hf_div(A, O), A));
    }
    // Multiplicative inverse of one: inv(1) = 1
    assert(hf_equal(hf_inv(hf_one()), hf_one()));
}

void test_double_negation() {
  for (uint8_t a = 0; a < 17; a++) {
    HF A = hf_new(a);
    HF negA = hf_neg(A);
    HF doubleNegA = hf_neg(negA);
    // Check -(-A) == A
    assert(hf_equal(A, doubleNegA));
  }
}

void test_f17_large_values() {
  // Test hf_new with values around multiples of the modulus
  assert(hf_equal(hf_new(17), hf_zero())); // 17 % 17 = 0
  assert(hf_equal(hf_new(-1), hf_new(16))); // -1 mod 17 = 16
  assert(hf_equal(hf_new(34), hf_new(0))); // 34 % 17 = 0
  assert(hf_equal(hf_new(51), hf_new(0))); // 51 % 17 = 0
  assert(hf_equal(hf_new(18), hf_new(1))); // 18 % 17 = 1
  assert(hf_equal(hf_new(-17), hf_zero())); // -17 mod 17 = 0
  assert(hf_equal(hf_new(-18), hf_new(16))); // -18 mod 17 -> -1 -> 16
  assert(hf_equal(hf_new(-35), hf_new(16))); // -35 mod 17 -> -1 -> 16

  // Test with values that might challenge int64_t if not handled carefully (though less relevant for mod 17)
  int64_t large_pos = 1700000000000000001; // = 17 * 1e17 + 1
  int64_t large_neg = -1700000000000000001; // = -(17 * 1e17 + 1)
  assert(hf_equal(hf_new(large_pos), hf_new(1)));
  assert(hf_equal(hf_new(large_neg), hf_new(16))); // -(1) mod 17 = 16

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
  test_field_axioms();
  test_identities_inverses();
  test_double_negation();
  test_f17_large_values();
  return 0;
}
