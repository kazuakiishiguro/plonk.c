#include <assert.h>
#include "g1.h"

#define ASSERT_G1(p1, p2)			      \
  do {                                                \
    assert((p1).x.value == (p2).x.value);	      \
    assert((p1).y.value == (p2).y.value);	      \
  } while (0)

void test_g1_vectors() {
  G1 g = g1_generator();
  G1 neg_g = g1_neg(&g);
  G1 two = g1_add(&g, &g);
  G1 neg_two = g1_neg(&two);
  G1 three = g1_add(&two, &g);
  G1 four = g1_add(&two, &two);
  G1 neg_four = g1_neg(&four);
  G1 five = g1_add(&four, &g);
  G1 six = g1_add(&five, &g);
  G1 eight = g1_add(&four, &four);
  G1 neg_eight = g1_neg(&eight);
  G1 nine = g1_add(&eight, &g);
  G1 sixteen = g1_add(&eight, &eight);
  G1 neg_sixteen = g1_neg(&sixteen);

  ASSERT_G1(g1_new(1, 99), neg_g);
  ASSERT_G1(g1_new(68, 74), two);
  ASSERT_G1(g1_new(68, 27), neg_two);
  ASSERT_G1(g1_new(26, 45), three);
  ASSERT_G1(g1_new(65, 98), four);
  ASSERT_G1(g1_new(65, 3), neg_four);
  ASSERT_G1(g1_new(12, 32), five);
  ASSERT_G1(g1_new(18, 49), eight);
  ASSERT_G1(g1_new(18, 52), neg_eight);
  ASSERT_G1(g1_new(18, 52), nine);
  ASSERT_G1(g1_new(1, 99), sixteen);
  ASSERT_G1(g1_new(1, 2), neg_sixteen);
  ASSERT_G1(g1_mul(&g, 1), g);
  ASSERT_G1(g1_mul(&g, 2), two);
  ASSERT_G1(g1_mul(&g, 2), two);
  ASSERT_G1(g1_mul(&g, 6), six);
}

int main() {
  test_g1_vectors();
  return 0;
}
