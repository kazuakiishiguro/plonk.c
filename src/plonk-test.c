#include <assert.h>
#include "plonk.h"

void test_plonk() {
  // step 1: create trustet setup (srs)
  // toxic waste
  GF secret = gf_new(2);
  // number of points in srs
  size_t n = 10;

  SRS srs = srs_create(secret, n);

  // step 2: create the plonk instance
  // number of omega powers (size of H)
  size_t num_constraints = 4;
  size_t omega_pows = num_constraints + 1;
  PLONK plonk = plonk_new(srs, omega_pows);

  // step 3: set up constraints and assignments
  // define the gates
  // gates:
  // - gate 1: mul a * b
  // - gate 2: mul a * b
  // - gate 3: mul a * b
  // - gate 4: mul a + b

  CONSTRAINTS constraints;
  constraints.num_constraints = num_constraints;

  // initialize q_m, q_l, q_r, q_o, q_c arrays
  constraints.q_m = (HF *)malloc(num_constraints * sizeof(HF));
  constraints.q_l = (HF *)malloc(num_constraints * sizeof(HF));
  constraints.q_r = (HF *)malloc(num_constraints * sizeof(HF));
  constraints.q_o = (HF *)malloc(num_constraints * sizeof(HF));
  constraints.q_c = (HF *)malloc(num_constraints * sizeof(HF));

  if (!constraints.q_m || !constraints.q_l || !constraints.q_r || !constraints.q_o || !constraints.q_c) {
    fprintf(stderr, "Memory allocation failed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // date definitions (for simplicity, we can define the gates manually)
  // for multiplication gates: q_m = l, q_o = -l
  // for addition gate: q_l = l, q_r = l, q_o = -l

  // gate 1: mul a * b
  constraints.q_m[0] = hf_one();
  constraints.q_l[0] = hf_zero();
  constraints.q_r[0] = hf_zero();
  constraints.q_o[0] = hf_neg(hf_one());
  constraints.q_c[0] = hf_zero();

  // gate 2: mul a * b
  constraints.q_m[1] = hf_one();
  constraints.q_l[1] = hf_zero();
  constraints.q_r[1] = hf_zero();
  constraints.q_o[1] = hf_neg(hf_one());
  constraints.q_c[1] = hf_zero();

  // gate 3: mul a * b
  constraints.q_m[2] = hf_one();
  constraints.q_l[2] = hf_zero();
  constraints.q_r[2] = hf_zero();
  constraints.q_o[2] = hf_neg(hf_one());
  constraints.q_c[2] = hf_zero();

  // gate 4: sum a + b
  constraints.q_m[3] = hf_zero();
  constraints.q_l[3] = hf_one();
  constraints.q_r[3] = hf_one();
  constraints.q_o[3] = hf_neg(hf_one());
  constraints.q_c[3] = hf_zero();

  // copy constraints (permutation)
  // for simplicity, let's assume c_a, c_b, c_c are arrays of COPY_OF structures
  constraints.c_a = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));
  constraints.c_b = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));
  constraints.c_c = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));

  if (!constraints.c_a || !constraints.c_b || !constraints.c_c) {
    fprintf(stderr, "Memory allocation failed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // define copy constraints
  // c_a: [b(1), b(2), b(3), c(1)]
  constraints.c_a[0] = (COPY_OF){COPYOF_B, 1};
  constraints.c_a[1] = (COPY_OF){COPYOF_B, 2};
  constraints.c_a[2] = (COPY_OF){COPYOF_B, 3};
  constraints.c_a[3] = (COPY_OF){COPYOF_C, 1};

  // c_b: [a(1), a(2), a(3), c(2)]
  constraints.c_b[0] = (COPY_OF){COPYOF_A, 1};
  constraints.c_b[1] = (COPY_OF){COPYOF_A, 2};
  constraints.c_b[2] = (COPY_OF){COPYOF_A, 3};
  constraints.c_b[3] = (COPY_OF){COPYOF_C, 2};

  // c_c: [a(4), a(4), a(4), c(3)]
  constraints.c_c[0] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[1] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[2] = (COPY_OF){COPYOF_A, 4};
  constraints.c_c[3] = (COPY_OF){COPYOF_C, 3};

  // step 4: define assignments
  ASSIGNMENTS assignments;
  assignments.len = num_constraints;
  assignments.a = (HF *)malloc(num_constraints * sizeof(HF));
  assignments.b = (HF *)malloc(num_constraints * sizeof(HF));
  assignments.c = (HF *)malloc(num_constraints * sizeof(HF));

  if (!assignments.a || !assignments.b || !assignments.c) {
    fprintf(stderr, "Memory allocation failed in test_plonk\n");
    exit(EXIT_FAILURE);
  }

  // assign values
  // assignment 1: a = 3, b = 3, c = 9
  assignments.a[0] = hf_new(3);
  assignments.b[0] = hf_new(3);
  assignments.c[0] = hf_new(9);

  // assignment 2: a = 4, b = 4, c = 16
  assignments.a[1] = hf_new(4);
  assignments.b[1] = hf_new(4);
  assignments.c[1] = hf_new(16);

  // assignment 3: a = 5, b = 5, c = 25
  assignments.a[2] = hf_new(5);
  assignments.b[2] = hf_new(5);
  assignments.c[2] = hf_new(25);

  // assignment 4: a = 9, b = 16, c = 25
  assignments.a[3] = hf_new(9);
  assignments.b[3] = hf_new(16);
  assignments.c[3] = hf_new(25);

  // step 5: define random numbers (the b's)
  HF rand[9] = {
    hf_new(7),
    hf_new(4),
    hf_new(11),
    hf_new(12),
    hf_new(16),
    hf_new(2),
    hf_new(14),
    hf_new(11),
    hf_new(7)
  };

  // step 6: define challenge calues
  CHALLENGE challenge;
  challenge.alpha = hf_new(15);
  challenge.beta = hf_new(12);
  challenge.gamma = hf_new(13);
  challenge.z = hf_new(5);
  challenge.z = hf_new(12);

  // step 7: call plonk_prove to generate a proof
  PROOF proof = plonk_prove(&plonk, &constraints, &assignments, &challenge, rand);

  // step 8: define the expected proof (for comparison)
  PROOF expected;
  expected.a_s = g1_new(hf_new(91).value, hf_new(66).value);
  expected.b_s = g1_new(hf_new(26).value, hf_new(45).value);
  expected.c_s = g1_new(hf_new(91).value, hf_new(35).value);
  expected.z_s = g1_new(hf_new(32).value, hf_new(59).value);
  expected.t_lo_s = g1_new(hf_new(12).value, hf_new(32).value);
  expected.t_mid_s = g1_new(hf_new(26).value, hf_new(45).value);
  expected.t_hi_s = g1_new(hf_new(91).value, hf_new(66).value);
  expected.w_z_s = g1_new(hf_new(91).value, hf_new(35).value);
  expected.w_z_omega_s = g1_new(hf_new(65).value, hf_new(98).value);
  expected.a_z = hf_new(15);
  expected.b_z = hf_new(13);
  expected.c_z = hf_new(5);
  expected.s_sigma_1_z = hf_one();
  expected.s_sigma_2_z = hf_new(12);
  expected.r_z = hf_new(15);
  expected.z_omega_z = hf_new(15);

  // step 9: compare the generated proof with the expected proof
  // TODO: check proof equality

  // step 10: verify the proof
  HF rand_verifier[1] = { hf_new(4) };
  assert(plonk_verify(&plonk, &constraints, &proof, &challenge, rand_verifier));

  free(constraints.q_m);
  free(constraints.q_l);
  free(constraints.q_r);
  free(constraints.q_o);
  free(constraints.q_c);
  free(constraints.c_a);
  free(constraints.c_b);
  free(constraints.c_c);
  free(assignments.a);
  free(assignments.b);
  free(assignments.c);
  plonk_free(&plonk);
}

int main() {
  test_plonk();
  return 0;
}
