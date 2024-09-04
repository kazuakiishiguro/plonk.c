#include <stdio.h>
#include <assert.h>
#include "g1.h"

#define ASSERT_G1(p1, p2)			      \
    do {					      \
        assert((p1).x.value == (p2).x.value);	      \
        assert((p1).y.value == (p2).y.value);	      \
    } while (0)

int main() {
     g1_p g = g1_p_generator();
     g1_p two = g1_p_add(g, g);
     g1_p three = g1_p_add(two, g);
     g1_p four = g1_p_add(two, two);
     g1_p five = g1_p_add(four, g);
     g1_p eight = g1_p_add(four, four);
     g1_p nine = g1_p_add(eight, g);
     g1_p sixteen = g1_p_add(eight, eight);
     ASSERT_G1(g1_p_new(68, 74), two);
     ASSERT_G1(g1_p_new(26, 45), three);
     ASSERT_G1(g1_p_new(65, 98), four);
     ASSERT_G1(g1_p_new(12, 32), five);
     ASSERT_G1(g1_p_new(18, 49), eight);
     ASSERT_G1(g1_p_new(18, 52), nine);
     ASSERT_G1(g1_p_new(1, 99), sixteen);
     return 0;
}
