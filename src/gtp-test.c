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
     gtp pow_6 = gtp_new(f101(42), f101(49));
     pow_6 = gtp_pow(&pow_6, 6);
     ASSERT_GTP(pow_6, gtp_new(f101(97), f101(89)));
     gtp base = gtp_new(f101(93), f101(76));
     gtp neg = gtp_sub_assign(&base);
     gtp pow_101 = gtp_pow(&base, 101);
     gtp pow_102 = gtp_pow(&base, 102);
     gtp neg_mul_base = gtp_mul(&neg, &base);
     ASSERT_GTP(pow_101, neg);
     ASSERT_GTP(pow_102, neg_mul_base);
     gtp pow_600 = gtp_new(f101(68), f101(47));
     pow_600 = gtp_pow(&pow_600, 600);
     ASSERT_GTP(pow_600, gtp_new(f101(97), f101(89)));
}

int main() {
     test_gtp_vectors();
     return 0;
}
