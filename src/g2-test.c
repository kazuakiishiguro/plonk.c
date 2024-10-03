#include <assert.h>
#include "g2.h"

#define ASSERT_G2(p1, p2)			      \
  do {                                                \
    assert((p1).x.value == (p2).x.value);	      \
    assert((p1).y.value == (p2).y.value);	      \
  } while (0)

void test_g2_vectors() {
  G2 g = g2_generator();
  G2 two = g2_add(&g, &g);
  G2 three = g2_add(&two, &g);
  G2 four = g2_add(&two, &two);
  G2 six = g2_add(&four, &two);

  ASSERT_G2(g2_new(90, 82), two);
  ASSERT_G2(g2_add(&three, &g), four);
  ASSERT_G2(g2_mul(g, 6), six);
}

int main () {
  test_g2_vectors();
  return 0;
}
