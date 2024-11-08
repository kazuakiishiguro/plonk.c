#include <assert.h>
#include "plonk.h"

#define ASSERT_POLY(p1, p2)                              \
  do {                                                   \
    assert(p1.len == p2.len);                            \
    for (size_t i = 0; i < p1.len; i++)                  \
      assert(hf_equal(p1.coeffs[i], p2.coeffs[i]));      \
  } while (0)

bool matrix_equal(const MATRIX *a, const MATRIX *b) {
  if (a->m != b->m || a->n != b->n)
    return false;
  size_t size = a->n * a->n;
  for (size_t i = 0; i < size; i++) {
    if (!hf_equal(a->v[i], b->v[i])) {
      return false;
    }
  }
  return true;
}

void test_plonk_init() {
  // create the trusted setup
  GF secret = f101(2); // the toxic waste
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);

  assert(plonk.h_len == h_len);

  // check roots of unity(H) generation
  for (uint8_t i = 0; i < h_len; i ++) {
    assert(hf_equal(plonk.h[i], hf_pow(hf_new(OMEGA_VALUE), i)));
  }

  HF vals[] = {f17(13), f17(13), f17(13), f17(13), f17(13), f17(16), f17(4), f17(1), f17(13), f17(4), f17(13), f17(4), f17(13), f17(1), f17(4), f17(16)};
  MATRIX expected = matrix_new(vals, h_len, h_len);
  assert(matrix_equal(&plonk.h_pows_inv, &expected));

  matrix_free(&expected);
  plonk_free(&plonk);
}

void test_interpolate_at_h() {
  GF secret = f101(2);
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);
  HF vec[] = {f17(3), f17(4), f17(0), f17(0)}; // adds paddings
  POLY result = interpolate_at_h(&plonk, vec, 4);
  // 6+x+4x^2+9x^3
  HF coeffs[] = {f17(6), f17(1), f17(4), f17(9)};
  POLY expected = poly_new(coeffs, 4);
  ASSERT_POLY(result, expected);

  plonk_free(&plonk);
}

void test_copy_constraints() {
  GF secret = f101(2);
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);

  CONSTRAINTS constraints;
  constraints.num_constraints = h_len;
  constraints.c_a = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_b = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_c = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  if (!constraints.c_a || !constraints.c_b || !constraints.c_c) {
    fprintf(stderr, "Memory allocation fialed in test_plonk\n");
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
  constraints.c_c[1] = (COPY_OF){COPYOF_B, 4};
  constraints.c_c[2] = (COPY_OF){COPYOF_C, 4};
  constraints.c_c[3] = (COPY_OF){COPYOF_C, 3};

  HF *sigma_c_a = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  HF *sigma_c_b = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  HF *sigma_c_c = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  copy_constraints_to_roots(&plonk, constraints.c_a, constraints.num_constraints, sigma_c_a);
  copy_constraints_to_roots(&plonk, constraints.c_b, constraints.num_constraints, sigma_c_b);
  copy_constraints_to_roots(&plonk, constraints.c_c, constraints.num_constraints, sigma_c_c);

  HF expected_c_a[] = {f17(2), f17(8), f17(15), f17(3)};
  HF expected_c_b[] = {f17(1), f17(4), f17(16), f17(12)};
  HF expected_c_c[] = {f17(13), f17(9), f17(5), f17(14)};
  for (uint8_t i = 0; i < constraints.num_constraints; i++) {
    assert(hf_equal(sigma_c_a[i], expected_c_a[i]));
    assert(hf_equal(sigma_c_b[i], expected_c_b[i]));
    assert(hf_equal(sigma_c_c[i], expected_c_c[i]));
  }

  free(constraints.c_a);
  free(constraints.c_b);
  free(constraints.c_c);
  free(sigma_c_a);
  free(sigma_c_b);
  free(sigma_c_c);
  plonk_free(&plonk);
}

void test_plonk_prove_verify() {
  // create the trusted setup
  GF secret = f101(2); // the toxic waste
  size_t n = 6;
  size_t h_len = 4;
  SRS srs = srs_create(secret, n);
  PLONK plonk = plonk_new(srs, h_len);

  // set up constraints and assignments
  // define the gates
  // gates:
  // - gate 1: mul a * b
  // - gate 2: mul a * b
  // - gate 3: mul a * b
  // - gate 4: mul a + b

  CONSTRAINTS constraints;
  constraints.num_constraints = h_len;

  // initialize q_m, q_l, q_r, q_o, q_c arrays
  constraints.q_m = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_l = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_r = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_o = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_c = (HF *)malloc(constraints.num_constraints * sizeof(HF));
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
  constraints.c_a = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_b = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_c = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));

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
  constraints.c_c[1] = (COPY_OF){COPYOF_B, 4};
  constraints.c_c[2] = (COPY_OF){COPYOF_C, 4};
  constraints.c_c[3] = (COPY_OF){COPYOF_C, 3};

  // step 4: define assignments
  ASSIGNMENTS assignments;
  assignments.len = constraints.num_constraints;
  assignments.a = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  assignments.b = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  assignments.c = (HF *)malloc(constraints.num_constraints * sizeof(HF));

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

  /* step 6: define challenge calues */
  CHALLENGE challenge;
  challenge.alpha = hf_new(15);
  challenge.beta = hf_new(12);
  challenge.gamma = hf_new(13);
  challenge.z = hf_new(5);
  challenge.v = hf_new(12);

  // step 7: call plonk_prove to generate a proof
  PROOF _ = plonk_prove(&plonk, &constraints, &assignments, &challenge, rand);

  free(constraints.q_m);
  free(constraints.q_l);
  free(constraints.q_r);
  free(constraints.q_o);
  free(constraints.q_c);
  free(constraints.c_a);
  free(constraints.c_b);
  free(constraints.c_c);
  plonk_free(&plonk);
}

int main() {
  test_plonk_init();
  test_interpolate_at_h();
  test_copy_constraints();
  test_plonk_prove_verify();
}
