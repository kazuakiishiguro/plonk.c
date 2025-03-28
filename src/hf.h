#ifndef HF_H
#define HF_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

#define MODULO_HF 17

/**
 * @bries Represents an element in the finite field GF(17).
 *        Value is always kep in the range [0, 16].
 */
typedef struct {
  uint8_t value;
} HF;

/**
 * @brief Creates a new field element from an integer.
 *        Result is reduced modulo MODULO_HF. Handles negative inputs correctly.
 * @param value The integer value.
 * @return HF The field element representing value mod MODULO_HF.
 */
static inline HF hf_new(int64_t value) {
  int64_t tmp = value % MODULO_HF;
  if (tmp < 0) {
    tmp += MODULO_HF;
  }

  HF hf;
  hf.value = (uint8_t)tmp;
  return hf;
}

/**
 * @brief Alias for hf_new. Creates a field element from an integer.
 * @param n The integer value.
 * @return HF The field element representing n mod MODULO_HF.
 */
static inline HF f17(int64_t n) {
  return hf_new(n);
}

/**
 * @brief Returns the zero element of the field (additive identity).
 * @return HF The zero element {0}.d
 */
static inline HF hf_zero() {
  HF r = {0};
  return r;
}

/**
 * @brief Returns the one element of the field (multiplicative identity).
 * @return HF The one element {1}.
 */
static inline HF hf_one() {
  HF r = {1};
  return r;
}

/**
 * @brief Checks if two field elements are equal.
 * @param a The first field element.
 * @param b The second field element.
 * @return bool True if a and b have the same value, false otherwise.
 */
static inline bool hf_equal(HF a, HF b) {
  return a.value == b.value;
}

/**
 * @brief Adds two field elements.
 * @param a The first field element.
 * @param b The second field element.
 * @return HF The sum (a + b) mod MODULO_HF.
 */
static inline HF hf_add(HF a, HF b) {
  uint8_t sum = a.value + b.value;
  if (sum >= MODULO_HF) sum -= MODULO_HF;
  HF r = { sum };
  return r;
}

/**
 * @brief Subtracts one field element from another.
 * @param a The minuend.
 * @param b The subtrahend.
 * @return HF The difference (a - b) mod MODULO_HF.
 */
static inline HF hf_sub(HF a, HF b) {
  int8_t diff = (int8_t)a.value - (int8_t)b.value;
  if (diff < 0) diff += MODULO_HF;
  HF r = { (uint8_t)diff };
  return r;
}

/**
 * @brief Multiplies two field elements.
 * @param a The first field element.
 * @param b The second field element.
 * @return HF The product (a * b) mod MODULO_HF.
 */
static inline HF hf_mul(HF a, HF b) {
  uint8_t product = (uint8_t)((uint16_t)a.value * (uint16_t)b.value % MODULO_HF);
  HF r = { product };
  return r;
}

/**
 * @brief Computes the additive inverse (negation) of a field element.
 * @param a The field element.
 * @return HF The element -a mod MODULO_HF.
 */
static inline HF hf_neg(HF a) {
  HF r = { a.value == 0 ? 0 : (uint8_t)(MODULO_HF - a.value) };
  return r;
}

/**
 * @brief Computes the power of a field element (exponentiation by squaring).
 * @param base The base field element.
 * @param exp The non-negative integer exponent.
 * @return HF The result (base ^ exp) mod MODULO_HF.
 */
static inline HF hf_pow(HF base, uint64_t exp) {
  HF r = hf_one();
  while (exp > 0) {
    if (exp & 1) {
      r = hf_mul(r, base);
    }
    base = hf_mul(base, base);
    exp >>= 1;
  }
  return r;
}

/**
 * @brief Lookup table for multiplicative inverses modulo MODULO_HF (17).
 *        `hf_inverses[i]` contains the value `j` such that `(i * j) % 17 == 1`.
 *        Index 0 is undefined/unused for inversion, typically stores 0.
 *        This turns hf_inv() into an O(1) lookup rather than doing exponentiation.
 */
static const uint8_t hf_inverses[17] = {
  // 0^-1 is undefined, store 0
  0,
  // 1*1=1 mod 17
  1,
  // 2*9=18=1 mod 17
  9,
  // 3*6=18=1 mod 17
  6,
  // 4*13=52=1 mod 17
  13,
  // 5*7=35=1 mod 17
  7,
  // 6*3=18=1 mod 17
  3,
  // 7*5=35=1 mod 17
  5,
  // 8*15=120=1 mod 17 (17*7 = 119)
  15,
  // 9*2=18=1 mod 17
  2,
  // 10*12=120=1 mod 17
  12,
  // 11*14=154=1 mod 17 (17*9 = 153)
  14,
  // 12*10=120=1 mod 17
  10,
  // 13*4=52=1 mod 17
  4,
  // 14*11=154=1 mod 17
  11,
  // 15*8=120=1 mod 17
  8,
  // 16*16 = (-1)*(-1) = 1 mod 17
  16
};

/**
 * @brief Computes the multiplicative inverse of a field element using a lookup table.
 * @param field The field element to invert. Assumed non-zero for mathematical correctness.
 *              If field is zero, returns 0 based on the table.
 * @return HF The multiplicative inverse (1 / field) mod MODULO_HF.
 */
static inline HF hf_inv(HF field) {
  HF r = { hf_inverses[field.value] };
  return r;
}

/**
 * @brief Divides one field element by another (a * b^-1).
 * @param a The numerator.
 * @param b The denominator. Must be non-zero for defined result.
 *          If b is zero, the result depends on hf_inv(0) which returns 0,
 *          leading to a result of a * 0 = 0.
 * @return HF The result (a / b) mod MODULO_HF.
 */
static inline HF hf_div(HF a, HF b) {
  return hf_mul(a, hf_inv(b));
}

#endif // HF_H
