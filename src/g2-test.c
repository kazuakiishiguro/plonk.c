#include <assert.h>
#include "g2.h"

#define ASSERT_G2(p1, p2)			      \
    do {					      \
        assert((p1).a.value == (p2).a.value);	      \
        assert((p1).b.value == (p2).b.value);	      \
    } while (0)

void test_g2_vectors() {
     g2_p g = g2_p_generator();
     g2_p two = g2_p_add(&g, &g);
     g2_p three = g2_p_add(&two, &g);
     g2_p four = g2_p_add(&two, &two);
     g2_p six = g2_p_add(&four, &two);

     ASSERT_G2(g2_p_new(90, 82), two);
     ASSERT_G2(g2_p_add(&three, &g), four);
     ASSERT_G2(g2_p_mul(g, f101(6)), six);
}

int main () {
     test_g2_vectors();
     return 0;
}
