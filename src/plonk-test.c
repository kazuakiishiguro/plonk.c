#include "plonk.h"

void test_plonk() {
  // step 1: create trustet setup (srs)
  // toxic waste
  u8_fe secret = u8_fe_new(2);
  // number of points in srs
  size_t n = 6;

  srs s = srs_create(secret, n);

  // step 2: create the plonk instance
  // number of omega powers (size of H)
  size_t omega_pows = 4;
  plonk p = plonk_new(s, omega_pows);

  // step 3: set up constraints and assignments
  // define the gates
  // gates:
  // - gate 1: mul a * b
  // - gate 2: mul a * b
  // - gate 3: mul a * b
  // - gate 4: mul a + b

  size_t num_constraints = 4;

  constraints c;
  c.num_constraints = num_constraints;

  // initialize q_m, q_l, q_r, q_o, q_c arrays
  c.q_m = (u8_fe *)malloc(num_constraints * sizeof(u8_fe));
  c.q_l = (u8_fe *)malloc(num_constraints * sizeof(u8_fe));
  c.q_r = (u8_fe *)malloc(num_constraints * sizeof(u8_fe));
  c.q_o = (u8_fe *)malloc(num_constraints * sizeof(u8_fe));
  c.q_c = (u8_fe *)malloc(num_constraints * sizeof(u8_fe));

  if (!c.q_m || !c.q_l || !c.q_r || !c.q_o || !c.q_c) {
    fprintf(stderr, "Memory allocation failed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // date definitions (for simplicity, we can define the gates manually)
  // for multiplication gates: q_m = l, q_o = -l
  // for addition gate: q_l = l, q_r = l, q_o = -l

  // gate 1: mul a * b
  c.q_m[0] = u8_fe_new(1);
  c.q_l[0] = u8_fe_new(0);
  c.q_r[0] = u8_fe_new(0);
  c.q_o[0] = u8_fe_neg(u8_fe_new(1));
  c.q_c[0] = u8_fe_new(0);

  // gate 2: mul a * b
  c.q_m[1] = u8_fe_new(1);
  c.q_l[1] = u8_fe_new(0);
  c.q_r[1] = u8_fe_new(0);
  c.q_o[1] = u8_fe_neg(u8_fe_new(1));
  c.q_c[1] = u8_fe_new(0);

  // gate 3: mul a * b
  c.q_m[2] = u8_fe_new(1);
  c.q_l[2] = u8_fe_new(0);
  c.q_r[2] = u8_fe_new(0);
  c.q_o[2] = u8_fe_neg(u8_fe_new(1));
  c.q_c[2] = u8_fe_new(0);

  // gate 4: sum a + b
  c.q_m[3] = u8_fe_new(0);
  c.q_l[3] = u8_fe_new(1);
  c.q_r[3] = u8_fe_new(1);
  c.q_o[3] = u8_fe_neg(u8_fe_new(1));
  c.q_c[3] = u8_fe_new(0);

  // copy constraints (permutation)
  // for simplicity, let's assume c_a, c_b, c_c are arrays of copy_of structures
  c.c_a = (copy_of *)malloc(num_constraints * sizeof(copy_of));
  c.c_b = (copy_of *)malloc(num_constraints * sizeof(copy_of));
  c.c_c = (copy_of *)malloc(num_constraints * sizeof(copy_of));

  if (!c.c_a || !c.c_b || !c.c_c) {
    fprintf(stderr, "Memory allocation failed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // define copy constraints
  // c_a: [b(1), b(2), b(3), c(1)]
  c.c_a[0] = (copy_of){COPYOF_B, 1};
  c.c_a[1] = (copy_of){COPYOF_B, 2};
  c.c_a[2] = (copy_of){COPYOF_B, 3};
  c.c_a[3] = (copy_of){COPYOF_C, 1};

  // c_b: [a(1), a(2), a(3), c(2)]
  c.c_b[0] = (copy_of){COPYOF_A, 1};
  c.c_b[1] = (copy_of){COPYOF_A, 2};
  c.c_b[2] = (copy_of){COPYOF_A, 3};
  c.c_b[3] = (copy_of){COPYOF_C, 2};

  // c_c: [a(4), a(4), a(4), c(3)]
  c.c_c[0] = (copy_of){COPYOF_A, 4};
  c.c_c[1] = (copy_of){COPYOF_A, 4};
  c.c_c[2] = (copy_of){COPYOF_A, 4};
  c.c_c[3] = (copy_of){COPYOF_C, 3};
}

int main() {
  return 0;
}
