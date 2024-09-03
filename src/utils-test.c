#include <assert.h>
#include <stdio.h>
#include "utils.h"

int main() {
    uint64_t value = 2;
    uint64_t modulo = 13;
    u64_fe field = u64_fe_new(value, modulo);
    u64_fe result = u64_fe_pow(&field, 5);
    assert(result.value == 6);

    u64_fe inv = u64_fe_inv(&field);
    assert(inv.value == 10);
    return 0;
}
