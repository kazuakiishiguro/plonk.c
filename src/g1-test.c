#include <assert.h>
#include "g1.h"

#define ASSERT_G1(p1, p2)			      \
    do {					      \
        assert((p1).x.value == (p2).x.value);	      \
        assert((p1).y.value == (p2).y.value);	      \
    } while (0)

void test_g1_vectors() {
     g1_p g = g1_p_generator();
     g1_p neg_g = g1_p_neg(&g);
     g1_p two = g1_p_add(&g, &g);
     g1_p neg_two = g1_p_neg(&two);
     g1_p three = g1_p_add(&two, &g);
     g1_p four = g1_p_add(&two, &two);
     g1_p neg_four = g1_p_neg(&four);
     g1_p five = g1_p_add(&four, &g);
     g1_p six = g1_p_add(&five, &g);
     g1_p eight = g1_p_add(&four, &four);
     g1_p neg_eight = g1_p_neg(&eight);
     g1_p nine = g1_p_add(&eight, &g);
     g1_p sixteen = g1_p_add(&eight, &eight);
     g1_p neg_sixteen = g1_p_neg(&sixteen);

     ASSERT_G1(g1_p_new(1, 99), neg_g);
     ASSERT_G1(g1_p_new(68, 74), two);
     ASSERT_G1(g1_p_new(68, 27), neg_two);
     ASSERT_G1(g1_p_new(26, 45), three);
     ASSERT_G1(g1_p_new(65, 98), four);
     ASSERT_G1(g1_p_new(65, 3), neg_four);
     ASSERT_G1(g1_p_new(12, 32), five);
     ASSERT_G1(g1_p_new(18, 49), eight);
     ASSERT_G1(g1_p_new(18, 52), neg_eight);
     ASSERT_G1(g1_p_new(18, 52), nine);
     ASSERT_G1(g1_p_new(1, 99), sixteen);
     ASSERT_G1(g1_p_new(1, 2), neg_sixteen);
     ASSERT_G1(g1_p_mul(g, f101(1)), g);
     ASSERT_G1(g1_p_mul(g, f101(2)), two);
     ASSERT_G1(g1_p_mul(g, f101(2)), two);
     ASSERT_G1(g1_p_mul(g, f101(6)), six);
}

int main() {
     test_g1_vectors();
     return 0;
}
