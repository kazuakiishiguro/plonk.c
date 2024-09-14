#include <assert.h>
#include "gtp.h"

#define ASSERT_GTP(p1, p2)			      \
    do {					      \
        assert((p1).a.value == (p2).a.value);	      \
        assert((p1).b.value == (p2).b.value);	      \
    } while (0)

void test_gtp_vectors() {
     gtp a = gtp_new(f101(26), f101(97));
     gtp b = gtp_new(f101(93), f101(76));
     ASSERT_GTP(gtp_mul(&a, &b), gtp_new(f101(97), f101(89)));
}

int main() {
     test_gtp_vectors();
     return 0;
}
