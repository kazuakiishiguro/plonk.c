#ifndef FE_H
#define FE_H

#include <stdbool.h>
#include <stdint.h>
#include "hf.h"

#define MODULO_GF 101

/**
 * @brief Represents an element in the finite field GF(101).
 *        Value is always kept in the range [0, 100].
 */
typedef struct {
  uint8_t value;
} GF;

/**
 * @brief Creates a new field element from an integer.
 *        Result is reduced modulo MODULO_GF. Handles negative inputs correctly.
 * @param value The integer value.
 * @return GF The field element representing value mod MODULO_GF.
 */
static inline GF gf_new(int64_t value) {
  int64_t tmp = value % MODULO_GF; // Note: % operator behavir for negative numbers depends on C standard/compiler
  if (tmp < 0) {
    tmp += MODULO_GF;
  }

  GF gf;
  gf.value = (uint8_t)tmp;
  return gf;
}

/**
 * @brief Alias for gf_new. Creates a field element from an integer.
 * @param n The integer value.
 * @return GF The field element representing n mod MODULO_GF.
 */
static inline GF f101(int64_t n) {
  return gf_new(n);
}

/**
 * @brief Returns the zero element of the field (additive identity).
 * @return GF The zero element {0}.
 */
static inline GF gf_zero() {
  GF r = {0};
  return r;
}

/**
 * @brief Checks if a number is odd. Helper function.
 * @param n The number to check.
 * @return bool True if n is odd, false otherwise.
 */
static inline bool is_odd(uint64_t n) {
  return n & 1;
}

/**
 * @brief Returns the one element of the field (multiplicative identity).
 * @return GF The one element {1}.
 */
static inline GF gf_one() {
  GF r = {1};
  return r;
}

/**
 * @brief Checks if two field elements are equal.
 * @param a The first field element.
 * @param b The second field element.
 * @return bool True if a and b have the same value, false otherwise.
 */
static inline bool gf_equal(GF a, GF b) {
  return a.value == b.value;
}

/**
 * @brief Adds two field elements.
 * @param a The first field element.
 * @param b The second field element.
 * @return GF The sum (a + b) mod MODULO_GF.
 */
static inline GF gf_add(GF a, GF b) {
  uint16_t sum = a.value + b.value;
  if (sum >= MODULO_GF) sum -= MODULO_GF;
  GF r;
  r.value = (uint8_t)sum; // This cast is safe because sum < 2*MODULO_GF
  return r;
}

/**
 * @brief Subtracts one field element from another.
 * @param a The minuend.
 * @param b The subtrahend.
 * @return GF The difference (a - b) mod MODULO_GF.
 */
static inline GF gf_sub(GF a, GF b) {
  int16_t diff = (int16_t)a.value - (int16_t)b.value;
  if (diff < 0) diff += MODULO_GF;
  GF r;
  r.value = (uint8_t)diff; // Safe cast as 0 <= diff < MODULO_GF
  return r;
}

/**
 * @brief Multiplies two field elements.
 * @param a The first field element.
 * @param b The second field element.
 * @return HF The product (a * b) mod MODULO_GF.
 */
static inline GF gf_mul(GF a, GF b) {
  uint16_t product = (uint16_t)a.value * (uint16_t)b.value;
  GF r;
  r.value = product % MODULO_GF;
  return r;
}

/**
 * @brief Computes the additive inverse (negation) of a field element.
 * @param a The field element.
 * @return GF The element -a mod MODULO_GF.
 */
static inline GF gf_neg(GF a) {
  if (a.value == 0) return gf_zero();
  GF r;
  r.value = a.value == 0 ? 0 : MODULO_GF - a.value;
  return r;
}

/**
 * @brief Computes the power of a field element (exponentiation by squaring).
 * @param base The base field element.
 * @param exp The non-negative integer exponent.
 * @return GF The result (base ^ exp) mod MODULO_GF.
 */
static inline GF gf_pow(GF field, uint64_t exp) {
  GF r = gf_one();
  GF base = field;
  while (exp > 0) {
    if (is_odd(exp)) {
      r = gf_mul(r, base);
    }
    exp >>= 1;
    base = gf_mul(base, base);
  }
  return r;
}

/**
 * @brief Computes the multiplicative inverse using Fermat's Little Theorem.
 *        Requires MODULO_GF to be prime. a^(p-2) mod p = a^-1 mod p.
 * @param field The field element to invert. Assumed non-zero. Behavior for 0 is 0^99 = 0.
 * @return GF The multiplicative inverse (1 / field) mod MODULO_GF.
 */
static inline GF gf_inv(GF field) {
  // Using Fermat's Little Theorem: a^(p-2) mod p
  return gf_pow(field, MODULO_GF -2);
}

/**
 * @brief Divides one field element by another (a * b^-1).
 * @param a The numerator.
 * @param b The denominator. Must be non-zero for defined result mathematically.
 * @return GF The result (a / b) mod MODULO_GF. If b is zero, result depends on gf_inv(0).
 */
static inline GF gf_div(GF a, GF b) {
  return gf_mul(a, gf_inv(b));
}

/**
 * @brief Converts an element from GF(17) (HF) to GF(101) (GF).
 *        This is a direct value conversion followed by reduction mod 101.
 * @param hf_elem The GF(17) element.
 * @return GF The corresponding element in GF(101).
 */
static inline GF gf_from_hf(HF hf_elem) {
  return gf_new(hf_elem.value);
}

#endif
