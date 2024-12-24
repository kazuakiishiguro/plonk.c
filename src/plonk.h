#ifndef PLONK_H
#define PLONK_H

#include <assert.h>
#include <stdlib.h>
#include "constraints.h"
#include "matrix.h"
#include "pairing.h"
#include "poly.h"
#include "srs.h"
#include "hf.h"

#define OMEGA_VALUE 4 // example value; should be a primitive roor of unity in the field
#define K1_VALUE    2 // should not be in H
#define K2_VALUE    3 // should not be in H or K1 * H

typedef struct {
  HF alpha;
  HF beta;
  HF gamma;
  HF z;
  HF v;
} CHALLENGE;

typedef struct {
  G1 a_s;
  G1 b_s;
  G1 c_s;
  G1 z_s;
  G1 t_lo_s;
  G1 t_mid_s;
  G1 t_hi_s;
  G1 w_z_s;
  G1 w_z_omega_s;
  HF a_z;
  HF b_z;
  HF c_z;
  HF s_sigma_1_z;
  HF s_sigma_2_z;
  HF r_z;
  HF z_omega_z;
} PROOF;

typedef struct {
  SRS srs;
  HF *h;        // array of HF elements
  MATRIX h_pows_inv;
  size_t h_len;
  HF *k1_h;     // array of k1 * h elements
  HF *k2_h;     // array of k2 * h elements
  POLY z_h_x;   // polynomial z_h_x
} PLONK;

PLONK plonk_new(SRS srs, size_t n) {
  PLONK plonk;
  plonk.srs = srs;
  plonk.h_len = n;

  // init
  HF OMEGA = hf_new(OMEGA_VALUE);
  HF K1 = hf_new(K1_VALUE);
  HF K2 = hf_new(K2_VALUE);

  // compute h = [omega^0, omega^1, ...., omega^n]
  plonk.h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.h) {
    fprintf(stderr, "Memory allocation failed in plnk_new\n");
    exit(EXIT_FAILURE);
  }
  for (uint8_t i = 0; i < plonk.h_len; i++) {
    plonk.h[i] = hf_pow(OMEGA, i);
  }

  // check K1 and K2 are not in H
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.h[i], K1) || hf_equal(plonk.h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }

  // compute h1_h: h1_h[i] = K1 * h[i]
  // compute k2_h: k2_h[i] = K2 * h[i]
  plonk.k1_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  plonk.k2_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.k1_h || !plonk.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  // k1_h is a coset of H
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k1_h[i] = hf_mul(plonk.h[i], K1);
  }
  // check K2 is chosen so that it is neither an element of H nor k1_h[i]
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.k1_h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }
  // k2_h is a coset of H
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k2_h[i] = hf_mul(plonk.h[i], K2);
  }

  // build h_pows matrix (vandermonde matrix) and compute its inverse
  MATRIX h_pows = matrix_zero(plonk.h_len, plonk.h_len);
  for (size_t c = 0; c < h_pows.n; c++) {
    for (size_t r = 0; r < h_pows.m; r++) {
      matrix_set(&h_pows, r, c, hf_pow(plonk.h[r], c));
    }
  }
  plonk.h_pows_inv = matrix_inv(&h_pows);
  matrix_free(&h_pows);

  // compute z_h_x = Poly::z(h)
  plonk.z_h_x = poly_z(plonk.h, plonk.h_len);

  return plonk;
}

void plonk_free(PLONK *plonk) {
  srs_free(&plonk->srs);
  matrix_free(&plonk->h_pows_inv);
  poly_free(&plonk->z_h_x);

  if (plonk->h) {
    free(plonk->h);
    plonk->h = NULL;
  }

  if (plonk->k1_h) {
    free(plonk->k1_h);
    plonk->k1_h = NULL;
  }

  if (plonk->k2_h) {
    free(plonk->k2_h);
    plonk->k2_h = NULL;
  }
}

void copy_constraints_to_roots(const PLONK *plonk, const COPY_OF *copy_of, size_t len, HF *sigma) {
  for (size_t i = 0; i < len; i++) {
    size_t idx = copy_of[i].index - 1;
    switch (copy_of[i].type) {
    case COPYOF_A:
      sigma[i] = plonk->h[idx];
      break;
    case COPYOF_B:
      sigma[i] = plonk->k1_h[idx];
      break;
    case COPYOF_C:
      sigma[i] = plonk->k2_h[idx];
      break;
    default:
      fprintf(stderr, "Invalid copy_of type\n");
      exit(EXIT_FAILURE);
    }
  }
}

POLY interpolate_at_h(const PLONK *plonk, const HF *values, size_t len) {
  // ensure len == plonk->h_len
  if (len != plonk->h_len) {
    fprintf(stderr, "Length mismatch in interpolate_at_h: len= %zu, plonk->h_len = %zu\n", len, plonk->h_len);
    exit(EXIT_FAILURE);
  }

  MATRIX vec = matrix_zero(len, 1);
  for (size_t i = 0; i < len; i++) {
    matrix_set(&vec, i, 0, values[i]);
  }

  // multiply h_pows_inv * vec
  MATRIX res = matrix_mul(&plonk->h_pows_inv, &vec);

  // convert result matrix to polynomial
  HF *coeffs = (HF *)malloc(res.m * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }

  POLY result = poly_new(coeffs, res.m);

  // clean up
  matrix_free(&vec);
  matrix_free(&res);
  free(coeffs);

  return result;
}

void poly_print(const POLY *p) {
  bool first = true;
  for (size_t i = 0; i < p->len; i++) {
    HF coeff = p->coeffs[i];
    if (!hf_equal(coeff, hf_new(0))) {
      if (!first)
        printf("+");
      if (i == 0) {
        printf("%u", coeff.value);
      } else if (i == 1){
        printf("%ux", coeff.value);
      } else {
        if (coeff.value != 1)
          printf("%u", coeff.value);
        printf("x^%zu", i);
      }
      first = false;
    }
  }
  if (first) {
    printf("0");
  }
  printf("\n");
}

PROOF plonk_prove(
    PLONK *plonk,
    CONSTRAINTS *constraints,
    ASSIGNMENTS *assignments,
    CHALLENGE *challenge,
    HF rand[9]
) {
  PROOF proof;
  memset(&proof, 0, sizeof(PROOF)); // Initialize proof to zero
  bool success = false;             // We'll set true once everything is OK

  // -------------------------------
  // 1. Basic checks
  // -------------------------------
  // If an assert fails, the program exits immediately.
  // No memory is allocated yet, so no leak possible.
  assert(constraints_satisfy(constraints, assignments));

  // extract challenges
  HF alpha = challenge->alpha;
  HF beta = challenge->beta;
  HF gamma = challenge->gamma;
  HF z = challenge->z;
  HF v = challenge->v;

  // constants
  size_t n = constraints->num_constraints;
  HF omega = hf_new(OMEGA_VALUE);
  HF k1 = hf_new(K1_VALUE);
  HF k2 = hf_new(K2_VALUE);

  // -------------------------------
  // 2. Allocate sigmas
  // -------------------------------
  HF *sigma_1 = NULL;
  HF *sigma_2 = NULL;
  HF *sigma_3 = NULL;

  sigma_1 = (HF *)malloc(n * sizeof(HF));
  sigma_2 = (HF *)malloc(n * sizeof(HF));
  sigma_3 = (HF *)malloc(n * sizeof(HF));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_prove\n");
    goto cleanup;
  }
  copy_constraints_to_roots(plonk, constraints->c_a, n, sigma_1);
  copy_constraints_to_roots(plonk, constraints->c_b, n, sigma_2);
  copy_constraints_to_roots(plonk, constraints->c_c, n, sigma_3);

  // -------------------------------
  // 3. All polynomials we allocate
  // -------------------------------
  // We'll declare them as "poly_zero()" so we can safely poly_free() them
  // in the cleanup block even if they were never overwritten.
  //
  //   (a,b,c)                  : assignments
  //   (o,m,l,r,c)              : gate constraints
  //   (sigma1, sigma2, sigma3) : copy constraints
  POLY f_a_x = poly_zero();
  POLY f_b_x = poly_zero();
  POLY f_c_x = poly_zero();
  POLY q_o_x = poly_zero();
  POLY q_m_x = poly_zero();
  POLY q_l_x = poly_zero();
  POLY q_r_x = poly_zero();
  POLY q_c_x = poly_zero();
  POLY s_sigma_1x = poly_zero();
  POLY s_sigma_2x = poly_zero();
  POLY s_sigma_3x = poly_zero();

  // For round 1
  POLY a_x = poly_zero();
  POLY b_x = poly_zero();
  POLY c_x = poly_zero();

  // Accumulator polynomial
  POLY acc_x = poly_zero();

  // Z polynomial
  POLY z_x = poly_zero();

  // Lagrange
  POLY l_1_x = poly_zero();

  // Our t_x polynomial, plus partials
  POLY t_x_numer = poly_zero();
  POLY t_x = poly_zero();
  POLY remainder = poly_zero();

  // Sliced polynomials
  POLY t_lo_x = poly_zero();
  POLY t_mid_x = poly_zero();
  POLY t_hi_x = poly_zero();

  // R polynomials
  POLY r_1_x = poly_zero();
  POLY r_2_x = poly_zero();
  POLY r_3_x = poly_zero();
  POLY r_4_x = poly_zero();
  POLY r_x = poly_zero();

  // For final step
  POLY w_z_x = poly_zero();
  POLY w_z_x_quo = poly_zero();
  POLY w_z_omega_x = poly_zero();

  POLY t_mid_x_z = poly_zero();
  POLY t_hi_x_z = poly_zero();

  f_a_x = interpolate_at_h(plonk, assignments->a, plonk->h_len);
  f_b_x = interpolate_at_h(plonk, assignments->b, plonk->h_len);
  f_c_x = interpolate_at_h(plonk, assignments->c, plonk->h_len);
  q_o_x = interpolate_at_h(plonk, constraints->q_o, plonk->h_len);
  q_m_x = interpolate_at_h(plonk, constraints->q_m, plonk->h_len);
  q_l_x = interpolate_at_h(plonk, constraints->q_l, plonk->h_len);
  q_r_x = interpolate_at_h(plonk, constraints->q_r, plonk->h_len);
  q_c_x = interpolate_at_h(plonk, constraints->q_c, plonk->h_len);
  s_sigma_1x = interpolate_at_h(plonk, sigma_1, plonk->h_len);
  s_sigma_2x = interpolate_at_h(plonk, sigma_2, plonk->h_len);
  s_sigma_3x = interpolate_at_h(plonk, sigma_3, plonk->h_len);

  // Once we have s_sigma_1x, etc., we can free sigma arrays
  free(sigma_1); sigma_1 = NULL;
  free(sigma_2); sigma_2 = NULL;
  free(sigma_3); sigma_3 = NULL;

  // -------------------------------
  // 4. Round 1: a_x, b_x, c_x
  // -------------------------------
  HF b1 = rand[0], b2 = rand[1], b3 = rand[2], b4 = rand[3], b5 = rand[4], b6 = rand[5];

  // a_x = (b2 + b1*x) * z_h_x + f_a_x
  HF a_blinding_coeffs[] = {b2, b1};
  POLY a_blind_poly  = poly_new(a_blinding_coeffs, 2);
  POLY a_x_blinded = poly_mul(&a_blind_poly, &plonk->z_h_x);
  a_x = poly_add(&a_x_blinded, &f_a_x);

  // b_x = (b4 + b3*x) * z_h_x + f_b_x
  HF b_blinding_coeffs[] = {b4, b3};
  POLY b_blind_poly  = poly_new(b_blinding_coeffs, 2);
  POLY b_x_blinded = poly_mul(&b_blind_poly, &plonk->z_h_x);
  b_x = poly_add(&b_x_blinded, &f_b_x);

  // c_x = (b6 + b5*x) * z_h_x + f_c_x
  HF c_blinding_coeffs[] = {b6, b5};
  POLY c_blind_poly  = poly_new(c_blinding_coeffs, 2);
  POLY c_x_blinded = poly_mul(&c_blind_poly, &plonk->z_h_x);
  c_x = poly_add(&c_x_blinded, &f_c_x);

  // Evaluate at s
  G1 a_s = srs_eval_at_s(&plonk->srs, &a_x);
  G1 b_s = srs_eval_at_s(&plonk->srs, &b_x);
  G1 c_s = srs_eval_at_s(&plonk->srs, &c_x);

  // Free only the temporary polynomials used in building a_x, b_x, c_x
  poly_free(&a_blind_poly);
  poly_free(&a_x_blinded);
  poly_free(&b_blind_poly);
  poly_free(&b_x_blinded);
  poly_free(&c_blind_poly);
  poly_free(&c_x_blinded);

  // Also we can now free f_a_x, f_b_x, f_c_x if we’re done with them:
  poly_free(&f_a_x);
  poly_free(&f_b_x);
  poly_free(&f_c_x);

  // -------------------------------
  // 5. Round 2: accumulator vector
  // -------------------------------
  // check https://vitalik.ca/general/2019/09/22/plonk.html

  // initialize the accumulator vector with 1
  //
  //   beta   blinding against prover addition manipulation
  //   gamma  blinding against prover multiplication manipulation
  //   k1,k2  coset independance between a,b,c
  HF *acc = (HF *)malloc(n * sizeof(HF));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    goto cleanup;
  }

  acc[0] = hf_one();
  for (size_t i = 1; i < n; i++) {
    // assignments at positoin i - 1
    HF as_a = assignments->a[i - 1];
    HF as_b = assignments->b[i - 1];
    HF as_c = assignments->c[i - 1];

    // omega_pow is the root of unity
    HF omega_pow = hf_pow(omega, i - 1);

    // denominator and numerator computations
    HF denom = hf_mul(
        hf_mul(
            hf_add(as_a, hf_add(hf_mul(beta, omega_pow), gamma)),
            hf_add(as_b, hf_add(hf_mul(beta, hf_mul(k1, omega_pow)), gamma))
	       ),
        hf_add(as_c, hf_add(hf_mul(beta, hf_mul(k2, omega_pow)), gamma))
		      );

    // s_sigma evaluations at omega_pow
    HF s_sigma_1_eval = poly_eval(&s_sigma_1x, omega_pow);
    HF s_sigma_2_eval = poly_eval(&s_sigma_2x, omega_pow);
    HF s_sigma_3_eval = poly_eval(&s_sigma_3x, omega_pow);

    HF numer = hf_mul(
        hf_mul(
            hf_add(as_a, hf_add(hf_mul(beta, s_sigma_1_eval), gamma)),
            hf_add(as_b, hf_add(hf_mul(beta, s_sigma_2_eval), gamma))
	       ),
        hf_add(as_c, hf_add(hf_mul(beta, s_sigma_3_eval), gamma))
		      );

    HF fraction = hf_div(denom, numer);
    acc[i] = hf_mul(acc[i - 1], fraction);
  }

  // interpolate accumualtor vector
  acc_x = interpolate_at_h(plonk, acc, plonk->h_len);;
  free(acc);
  acc = NULL;

  // check that acc_x(omega^n) == 1
  HF omega_n = hf_pow(omega, n);
  HF acc_x_eval = poly_eval(&acc_x, omega_n);
  assert(hf_equal(acc_x_eval, hf_one()));

  // create z_x polynomial
  HF b7 = rand[6], b8 = rand[7], b9 = rand[8];
  HF z_blinding_coeffs[] = {b9, b8, b7};
  POLY z_blinding_poly = poly_new(z_blinding_coeffs, 3);
  POLY z_blinded = poly_mul(&z_blinding_poly, &plonk->z_h_x);
  z_x = poly_add(&z_blinded, &acc_x);

  // output of second step
  // evaluate z_x at s
  G1 z_s = srs_eval_at_s(&plonk->srs, &z_x);

  // clean up
  poly_free(&z_blinding_poly);
  poly_free(&z_blinded);


  // -------------------------------
  // 6. Compute t(x)
  // -------------------------------
  // build L1(x)
  HF *lagrange_vector = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed for lagrange_vector\n");
    exit(EXIT_FAILURE);
  }
  lagrange_vector[0] = hf_one(); // L1(omega^0) = 1, rest are 0
  l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);
  free(lagrange_vector);
  lagrange_vector = NULL;

  // compute p_i_x (public input polynomial)
  // assuming no public inputs for simplicity
  POLY p_i_x = poly_zero();

  // compute helper polynomials
  // compute t_1_z_h
  POLY a_x_b_x = poly_mul(&a_x, &b_x);
  POLY a_x_b_x_q_m_x = poly_mul(&a_x_b_x, &q_m_x);
  poly_free(&a_x_b_x);

  POLY a_x_q_l_x = poly_mul(&a_x, &q_l_x);
  POLY b_x_q_r_x = poly_mul(&b_x, &q_r_x);
  POLY c_x_q_o_x = poly_mul(&c_x, &q_o_x);

  POLY sum1 = poly_add(&a_x_b_x_q_m_x, &a_x_q_l_x);
  poly_free(&a_x_b_x_q_m_x);
  poly_free(&a_x_q_l_x);

  POLY sum2 = poly_add(&b_x_q_r_x, &c_x_q_o_x);
  poly_free(&b_x_q_r_x);
  poly_free(&c_x_q_o_x);

  POLY t_1_z_h = poly_add(&sum1, &sum2);
  poly_free(&sum1);
  poly_free(&sum2);
  t_1_z_h = poly_add(&t_1_z_h, &p_i_x);
  poly_free(&p_i_x);
  t_1_z_h = poly_add(&t_1_z_h, &q_c_x);

  // compute t_2_z_h
  HF beta_x_gamma_vec[] = {gamma, beta};
  POLY beta_x_gamma = poly_new(beta_x_gamma_vec, 2);
  POLY a_x_beta_x_gamma = poly_add(&a_x, &beta_x_gamma);
  POLY alpha_a_x_beta_x_gamma = poly_scale(&a_x_beta_x_gamma, alpha);
  HF beta_k1_x_gamma_vec[] = {gamma, hf_mul(beta, k1)};
  POLY beta_k1_x_gamma = poly_new(beta_k1_x_gamma_vec, 2);
  POLY b_x_beta_k1_x_gamma = poly_add(&b_x, &beta_k1_x_gamma);
  HF beta_k2_x_gamma_vec[] = {gamma, hf_mul(beta, k2)};
  POLY beta_k2_x_gamma = poly_new(beta_k2_x_gamma_vec, 2);
  POLY c_x_beta_k2_x_gamma = poly_add(&c_x, &beta_k2_x_gamma);
  POLY t_2_z_h = poly_mul(&alpha_a_x_beta_x_gamma, &b_x_beta_k1_x_gamma);
  t_2_z_h = poly_mul(&t_2_z_h, &c_x_beta_k2_x_gamma);
  t_2_z_h = poly_mul(&t_2_z_h, &z_x);

  poly_free(&beta_x_gamma);
  poly_free(&a_x_beta_x_gamma);
  poly_free(&beta_k1_x_gamma);
  poly_free(&beta_k2_x_gamma);
  poly_free(&alpha_a_x_beta_x_gamma);
  poly_free(&b_x_beta_k1_x_gamma);
  poly_free(&c_x_beta_k2_x_gamma);

  // compute t_3_z_h
  // alpha_a_x_beta_s_sigma1_x_gamma *
  // b_x_beta_s_sigma2_x_gamma *
  // c_x_teta_s_sigma3_x_gamma *
  // &z_omega_x
  POLY beta_s_sigma1 = poly_scale(&s_sigma_1x, beta);
  POLY a_x_beta_s_sigma1 = poly_add(&a_x, &beta_s_sigma1);
  POLY a_x_beta_s_sigma1_gamma = poly_add_hf(&a_x_beta_s_sigma1, gamma);
  POLY alpha_a_x_beta_s_sigma1_gamma = poly_scale(&a_x_beta_s_sigma1_gamma, alpha);
  POLY beta_s_sigma2 = poly_scale(&s_sigma_2x, beta);
  POLY b_x_beta_s_sigma2 = poly_add(&b_x, &beta_s_sigma2);
  POLY b_x_beta_s_sigma2_gamma = poly_add_hf(&b_x_beta_s_sigma2, gamma);

  POLY beta_s_sigma3 = poly_scale(&s_sigma_3x, beta);
  POLY c_x_beta_s_sigma3 = poly_add(&c_x, &beta_s_sigma3);
  POLY c_x_beta_s_sigma3_gamma = poly_add_hf(&c_x_beta_s_sigma3, gamma);

  // convert result matrix to polynomial
  HF *coeffs = (HF *)malloc(z_x.len * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    goto cleanup;
  }
  for (size_t i = 0; i < z_x.len; i++) {
    coeffs[i] = hf_mul(z_x.coeffs[i], hf_pow(omega, i));
  }
  POLY z_omega_x = poly_new(coeffs, z_x.len);
  free(coeffs);
  coeffs = NULL;

  POLY t_3_z_h = poly_mul(&alpha_a_x_beta_s_sigma1_gamma, &b_x_beta_s_sigma2_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &c_x_beta_s_sigma3_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &z_omega_x);

  poly_free(&beta_s_sigma1);
  poly_free(&a_x_beta_s_sigma1);
  poly_free(&beta_s_sigma2);
  poly_free(&b_x_beta_s_sigma2);
  poly_free(&beta_s_sigma3);
  poly_free(&c_x_beta_s_sigma3);

  // compute t_4_z_h
  HF neg_one[1];
  neg_one[0] = hf_neg(hf_one());
  POLY _1 = poly_new(neg_one, 1);
  POLY z_x_1 = poly_add(&z_x, &_1);
  POLY alhpa_2_z_x_1 = poly_scale(&z_x_1, hf_pow(alpha, 2));
  POLY t_4_z_h = poly_mul(&alhpa_2_z_x_1, &l_1_x);

  // build t_x_numer = t_1_z_h + t_2_z_h - t_3_z_h + t_4_z_h
  t_x_numer = poly_add(&t_1_z_h, &t_2_z_h);
  t_x_numer = poly_sub(&t_x_numer, &t_3_z_h);
  t_x_numer = poly_add(&t_x_numer, &t_4_z_h);

  poly_free(&t_1_z_h);
  poly_free(&t_2_z_h);
  poly_free(&t_3_z_h);
  poly_free(&t_4_z_h);
  poly_free(&_1);
  poly_free(&z_x_1);
  poly_free(&alhpa_2_z_x_1);

  poly_divide(&t_x_numer, &plonk->z_h_x, &t_x, &remainder);
  if (!poly_is_zero(&remainder)) {
    fprintf(stderr, "Non-zero remainder in t(x) division\n");
    poly_free(&remainder);
    goto cleanup;
  }
  poly_free(&remainder);

  // step 7: split t(x) into t_lo_x, t_mid_x, t_hi_x
  size_t degree = t_x.len;
  size_t part_size = n + 2;

  t_lo_x = poly_slice(&t_x, 0, part_size);
  t_mid_x = poly_slice(&t_x, part_size, 2 * part_size);
  t_hi_x = poly_slice(&t_x, 2 * part_size, degree);

  // step 8: evaluate t_lo_x, t_mid_x, t_hi_x at s
  G1 t_lo_s = srs_eval_at_s(&plonk->srs, &t_lo_x);
  G1 t_mid_s = srs_eval_at_s(&plonk->srs, &t_mid_x);
  G1 t_hi_s = srs_eval_at_s(&plonk->srs, &t_hi_x);

  // step 9: compute r(x) and evaluate it at z
  HF a_z = poly_eval(&a_x, z);
  HF b_z = poly_eval(&b_x, z);
  HF c_z = poly_eval(&c_x, z);
  HF s_sigma_1_z = poly_eval(&s_sigma_1x, z);
  HF s_sigma_2_z = poly_eval(&s_sigma_2x, z);
  HF t_z = poly_eval(&t_x, z);
  HF z_omega_z = poly_eval(&z_omega_x, z);
  poly_free(&t_x);

  // compute r_1_x
  POLY a_z_b_z_q_m_x = poly_scale(&q_m_x, hf_mul(a_z, b_z));
  POLY a_z_q_l_x = poly_scale(&q_l_x, a_z);
  POLY b_z_q_r_x = poly_scale(&q_r_x, b_z);
  POLY c_z_q_o_x = poly_scale(&q_o_x, c_z);
  r_1_x = poly_add(&a_z_b_z_q_m_x, &a_z_q_l_x);
  r_1_x = poly_add(&r_1_x, &b_z_q_r_x);
  r_1_x = poly_add(&r_1_x, &c_z_q_o_x);

  // compute r_2_x
  HF a_z_beta_z_gamma = hf_add(hf_add(a_z, hf_mul(beta, z)), gamma);
  HF b_z_beta_k1_z_gamma = hf_add(hf_add(b_z, hf_mul(hf_mul(beta, k1), z)), gamma);
  HF c_z_beta_k2_z_gamma = hf_add(hf_add(c_z, hf_mul(hf_mul(beta, k2), z)), gamma);
  r_2_x = poly_scale(&z_x,
                     hf_mul(
                         hf_mul(
                             hf_mul(a_z_beta_z_gamma, b_z_beta_k1_z_gamma),
                             c_z_beta_k2_z_gamma),
                         alpha));

  // compute r_3_x
  POLY s_sigma_3_beta_z_omega_z = poly_scale(&s_sigma_3x, hf_mul(beta, z_omega_z));
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma));
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma));
  r_3_x = poly_mul(&z_x, &s_sigma_3_beta_z_omega_z);
  r_3_x = poly_scale(&r_3_x,
                     hf_mul(
                         hf_mul(a_z_beta_s_sigma_1_z_gamma, b_z_beta_s_sigma_2_z_gamma),
                         alpha));

  // compute r_4_x
  r_4_x = poly_scale(&z_x, hf_mul(poly_eval(&l_1_x, z), hf_pow(alpha, 2)));

  r_x = poly_add(&r_1_x, &r_2_x);
  r_x = poly_add(&r_x, &r_3_x);
  r_x = poly_add(&r_x, &r_4_x);

  // compute linearization evaluation
  HF r_z = poly_eval(&r_x, z);

  poly_free(&q_o_x);
  poly_free(&q_m_x);
  poly_free(&q_l_x);
  poly_free(&q_r_x);
  poly_free(&q_c_x);

  // round 5
  // ---------------------------------------------------------------------------
  // we create two large polynomials that combine all the polynomials we've been
  // using so far and we output commitments to them.

  // compute opening proof polynomial w_z_x
  t_mid_x_z = poly_scale(&t_mid_x, hf_pow(z, n+2));
  t_hi_x_z = poly_scale(&t_hi_x, hf_pow(z, 2*n+4));
  w_z_x = poly_add(&t_lo_x, &t_mid_x_z);
  w_z_x = poly_add(&w_z_x, &t_hi_x_z);
  w_z_x = poly_add_hf(&w_z_x, hf_neg(t_z));
  POLY r_x_r_z_v = poly_add_hf(&r_x, hf_neg(r_z));
  r_x_r_z_v = poly_scale(&r_x_r_z_v, v);
  POLY a_x_a_z_v = poly_add_hf(&a_x, hf_neg(a_z));
  a_x_a_z_v = poly_scale(&a_x_a_z_v, hf_pow(v, 2));
  POLY b_x_b_z_v = poly_add_hf(&b_x, hf_neg(b_z));
  b_x_b_z_v = poly_scale(&b_x_b_z_v, hf_pow(v, 3));
  POLY c_x_c_z_v = poly_add_hf(&c_x, hf_neg(c_z));
  c_x_c_z_v = poly_scale(&c_x_c_z_v, hf_pow(v, 4));
  POLY s_sigma_1_1_z_v = poly_add_hf(&s_sigma_1x, hf_neg(s_sigma_1_z));
  s_sigma_1_1_z_v = poly_scale(&s_sigma_1_1_z_v, hf_pow(v, 5));
  POLY s_sigma_2_2_z_v = poly_add_hf(&s_sigma_2x, hf_neg(s_sigma_2_z));
  s_sigma_2_2_z_v = poly_scale(&s_sigma_2_2_z_v, hf_pow(v, 6));
  w_z_x = poly_add(&w_z_x, &r_x_r_z_v);
  w_z_x = poly_add(&w_z_x, &a_x_a_z_v);
  w_z_x = poly_add(&w_z_x, &b_x_b_z_v);
  w_z_x = poly_add(&w_z_x, &c_x_c_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_1_1_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_2_2_z_v);

  HF coeffs_denom1[] = {hf_neg(z), hf_one()};
  POLY denom1 = poly_new(coeffs_denom1, 2);
  POLY rem1;
  poly_divide(&w_z_x, &denom1, &w_z_x_quo, &rem1);
  assert(poly_is_zero(&rem1));

  POLY z_x_z_omega_z = poly_add_hf(&z_x, hf_neg(z_omega_z));
  HF coeffs_denom2[] = {hf_mul(hf_neg(z), omega), hf_one()};
  POLY denom2 = poly_new(coeffs_denom2, 2);
  POLY rem2;
  poly_divide(&z_x_z_omega_z, &denom2, &w_z_omega_x, &rem2);
  assert(poly_is_zero(&rem2));

  // compute opening proof polinomials at s
  G1 w_z_s = srs_eval_at_s(&plonk->srs, &w_z_x_quo);
  G1 w_z_omega_s = srs_eval_at_s(&plonk->srs, &w_z_omega_x);

  proof.a_s = a_s;
  proof.b_s = b_s;
  proof.c_s = c_s;
  proof.z_s = z_s;
  proof.t_lo_s = t_lo_s;
  proof.t_mid_s = t_mid_s;
  proof.t_hi_s = t_hi_s;
  proof.w_z_s = w_z_s;
  proof.w_z_omega_s = w_z_omega_s;
  proof.a_z = a_z;
  proof.b_z = b_z;
  proof.c_z = c_z;
  proof.s_sigma_1_z = s_sigma_1_z;
  proof.s_sigma_2_z = s_sigma_2_z;
  proof.r_z = r_z;
  proof.z_omega_z = z_omega_z;

  success = true;

cleanup:
  // -------------------------------
  // Free arrays
  // -------------------------------
  if (sigma_1) { free(sigma_1); sigma_1 = NULL; }
  if (sigma_2) { free(sigma_2); sigma_2 = NULL; }
  if (sigma_3) { free(sigma_3); sigma_3 = NULL; }

  poly_free(&t_lo_x);
  poly_free(&t_mid_x);
  poly_free(&t_hi_x);
  poly_free(&z_x);
  poly_free(&l_1_x);
  poly_free(&z_omega_x);
  poly_free(&s_sigma_1x);
  poly_free(&s_sigma_2x);
  poly_free(&s_sigma_3x);
  poly_free(&t_x_numer);
  poly_free(&r_1_x);
  poly_free(&r_2_x);
  poly_free(&r_3_x);
  poly_free(&r_4_x);
  poly_free(&r_x);
  poly_free(&t_mid_x_z);
  poly_free(&t_hi_x_z);

  // If not successful, zero out proof:
  if (!success) {
    memset(&proof, 0, sizeof(PROOF));
  }

  return proof;
}

bool plonk_verify(
    PLONK *plonk,
    CONSTRAINTS *constraints,
    PROOF *proof,
    CHALLENGE *challenge,
    HF rand[1]
) {
  bool success = false;
  size_t n = constraints->num_constraints;

  // ----------------------------------------------------
  // 0. Setup pointers/objects we'll allocate & need to free
  // ----------------------------------------------------
  HF *sigma_1 = NULL;
  HF *sigma_2 = NULL;
  HF *sigma_3 = NULL;

  // We'll store the polynomials we create in step 8 in these variables
  // so we can free them in case of an error.
  POLY q_m_x = poly_zero();
  POLY q_l_x = poly_zero();
  POLY q_r_x = poly_zero();
  POLY q_o_x = poly_zero();
  POLY q_c_x = poly_zero();
  POLY s_sigma_1_x = poly_zero();
  POLY s_sigma_2_x = poly_zero();
  POLY s_sigma_3_x = poly_zero();

  // We'll also create an extra "cleanup" label at the end
  // to free memory in all cases.

  // ----------------------------------------------------
  // 1. Unpack from proof and challenge
  // ---------------------------------------------------
  G1 a_s = proof->a_s;
  G1 b_s = proof->b_s;
  G1 c_s = proof->c_s;
  G1 z_s = proof->z_s;
  G1 t_lo_s = proof->t_lo_s;
  G1 t_mid_s = proof->t_mid_s;
  G1 t_hi_s = proof->t_hi_s;
  G1 w_z_s = proof->w_z_s;
  G1 w_z_omega_s = proof->w_z_omega_s;

  HF a_z = proof->a_z;
  HF b_z = proof->b_z;
  HF c_z = proof->c_z;
  HF s_sigma_1_z = proof->s_sigma_1_z;
  HF s_sigma_2_z = proof->s_sigma_2_z;
  HF r_z = proof->r_z;
  HF z_omega_z = proof->z_omega_z;

  // unpack challenge
  HF alpha = challenge->alpha;
  HF beta = challenge->beta;
  HF gamma = challenge->gamma;
  HF z = challenge->z;
  HF v = challenge->v;

  // constants
  HF omega = hf_new(OMEGA_VALUE);
  HF k1 = hf_new(K1_VALUE);
  HF k2 = hf_new(K2_VALUE);

  // ----------------------------------------------------
  // 2. Validate proof points in G1
  // ----------------------------------------------------
  if (!g1_is_on_curve(&a_s) ||
      !g1_is_on_curve(&b_s) ||
      !g1_is_on_curve(&c_s) ||
      !g1_is_on_curve(&z_s) ||
      !g1_is_on_curve(&t_lo_s) ||
      !g1_is_on_curve(&t_mid_s) ||
      !g1_is_on_curve(&t_hi_s) ||
      !g1_is_on_curve(&w_z_s) ||
      !g1_is_on_curve(&w_z_omega_s)) {
    // Return false immediately, no memory has been allocated yet, so no leak
    return false;
  }

  // ----------------------------------------------------
  // 3. Validate proof fields in HF
  // ----------------------------------------------------
  if (!hf_in_field(a_z) ||
      !hf_in_field(b_z) ||
      !hf_in_field(c_z) ||
      !hf_in_field(s_sigma_1_z) ||
      !hf_in_field(s_sigma_2_z) ||
      !hf_in_field(r_z) ||
      !hf_in_field(z_omega_z)) {
    // Also safe to return false now
    return false;
  }

  // ----------------------------------------------------
  // 4. Evaluate z_h at z
  // ----------------------------------------------------
  HF z_h_z = poly_eval(&plonk->z_h_x, z);

  // ----------------------------------------------------
  // 5. Evaluate Lagrange polynomial L1 at z
  // ----------------------------------------------------
  HF *lagrange_vector = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    goto cleanup;  // jump to free resources
  }

  lagrange_vector[0] = hf_new(1); // L1(omega^0) = 1
  POLY l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);

  // we can free lagrange_vector right away
  free(lagrange_vector);

  HF l_1_z = poly_eval(&l_1_x, z);
  poly_free(&l_1_x);

  // ----------------------------------------------------
  // 6. No public inputs => p_i_z = 0
  // ----------------------------------------------------
  HF p_i_z = hf_zero();

  // ----------------------------------------------------
  // 7. Compute quotient polynomial evaluation t_z
  // ----------------------------------------------------
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_1_z), gamma), a_z);
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_2_z), gamma), b_z);
  HF c_z_gamma = hf_add(c_z, gamma);
  HF l_1_z_alpha_2 = hf_mul(l_1_z, hf_pow(alpha, 2));

  // numerator = (r_z + p_i_z) - ( [a_z+beta*s_sigma_1_z+gamma]*[b_z+beta*s_sigma_2_z+gamma]*[c_z+gamma]*z_omega_z + alpha^2*l_1_z )
  HF numerator = hf_sub(hf_add(r_z, p_i_z),
                        hf_add(
                            hf_mul(
                                hf_mul(
                                    hf_mul(
                                        a_z_beta_s_sigma_1_z_gamma,
                                        b_z_beta_s_sigma_2_z_gamma),
                                    c_z_gamma),
                                z_omega_z),
                            l_1_z_alpha_2));

  // t_z = numerator / z_h_z
  HF t_z = hf_div(numerator, z_h_z);

  // ----------------------------------------------------
  // 8. Interpolate constraint polynomials & sigma polynomials
  //    Then evaluate them at s
  // ----------------------------------------------------
  // Allocate sigma arrays
  sigma_1 = (HF *)malloc(n * sizeof(HF));
  sigma_2 = (HF *)malloc(n * sizeof(HF));
  sigma_3 = (HF *)malloc(n * sizeof(HF));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    goto cleanup;
  }

  copy_constraints_to_roots(plonk, constraints->c_a, n, sigma_1);
  copy_constraints_to_roots(plonk, constraints->c_b, n, sigma_2);
  copy_constraints_to_roots(plonk, constraints->c_c, n, sigma_3);

  // Now build polynomials for q_m, q_l, q_r, q_o, q_c, s_sigma_1, s_sigma_2, s_sigma_3
  q_m_x = interpolate_at_h(plonk, constraints->q_m, n);
  q_l_x = interpolate_at_h(plonk, constraints->q_l, n);
  q_r_x = interpolate_at_h(plonk, constraints->q_r, n);
  q_o_x = interpolate_at_h(plonk, constraints->q_o, n);
  q_c_x = interpolate_at_h(plonk, constraints->q_c, n);
  s_sigma_1_x = interpolate_at_h(plonk, sigma_1, n);
  s_sigma_2_x = interpolate_at_h(plonk, sigma_2, n);
  s_sigma_3_x = interpolate_at_h(plonk, sigma_3, n);

  // Evaluate them at s
  G1 q_m_s = srs_eval_at_s(&plonk->srs, &q_m_x);
  G1 q_l_s = srs_eval_at_s(&plonk->srs, &q_l_x);
  G1 q_r_s = srs_eval_at_s(&plonk->srs, &q_r_x);
  G1 q_o_s = srs_eval_at_s(&plonk->srs, &q_o_x);
  G1 q_c_s = srs_eval_at_s(&plonk->srs, &q_c_x);
  G1 sigma_1_s = srs_eval_at_s(&plonk->srs, &s_sigma_1_x);
  G1 sigma_2_s = srs_eval_at_s(&plonk->srs, &s_sigma_2_x);
  G1 sigma_3_s = srs_eval_at_s(&plonk->srs, &s_sigma_3_x);

  // No further usage for sigma_1/2/3 arrays
  free(sigma_1);
  free(sigma_2);
  free(sigma_3);

  // ----------------------------------------------------
  // 9. Build partial commitments (d_1_s, d_2_s, d_3_s, etc.)
  // ----------------------------------------------------
  HF u = rand[0];

  // d_1_s = q_m_s*(a_z*b_z*v) + q_l_s*(a_z*v) + q_r_s*(b_z*v) + q_o_s*(c_z*v) + q_c_s*(v)
  // We'll do it step by step for clarity
  HF a_z_b_z_v = hf_mul(hf_mul(a_z, b_z), v);
  HF a_z_v = hf_mul(a_z, v);
  HF b_z_v = hf_mul(b_z, v);
  HF c_z_v = hf_mul(c_z, v);

  G1 tmp_m = g1_mul(&q_m_s, a_z_b_z_v.value);
  G1 tmp_l = g1_mul(&q_l_s, a_z_v.value);
  G1 tmp_r = g1_mul(&q_r_s, b_z_v.value);
  G1 tmp_o = g1_mul(&q_o_s, c_z_v.value);
  G1 tmp_c = g1_mul(&q_c_s, v.value);

  // sum them
  G1 d_1_s = g1_add(&tmp_m, &tmp_l);
  d_1_s = g1_add(&d_1_s, &tmp_r);
  d_1_s = g1_add(&d_1_s, &tmp_o);
  d_1_s = g1_add(&d_1_s, &tmp_c);

  // d_2_s = z_s * ( [a_z+beta*z+gamma]*[b_z+beta*k1*z+gamma]*[c_z+beta*k2*z+gamma]*alpha*v + l_1_z*alpha^2*v + u )
  HF term_d2_inner = hf_mul(
      hf_mul(
          hf_mul(
              hf_add(a_z, hf_add(hf_mul(beta, z), gamma)),
              hf_add(b_z, hf_add(hf_mul(beta, hf_mul(k1, z)), gamma))
                 ),
          hf_add(c_z, hf_add(hf_mul(beta, hf_mul(k2, z)), gamma))
             ),
      hf_mul(alpha, v)
                            );
  HF term_d2 = hf_add(term_d2_inner, hf_mul(l_1_z, hf_mul(hf_pow(alpha, 2), v)));
  term_d2 = hf_add(term_d2, u);
  G1 d_2_s = g1_mul(&z_s, term_d2.value);

  // d_3_s = sigma_3_s * ( [a_z + beta*s_sigma_1_z + gamma]*[b_z + beta*s_sigma_2_z + gamma]*alpha*v*beta*z_omega_z )
  HF term_d3_inner = hf_mul(
      hf_mul(
          hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma)),
          hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma))
             ),
      hf_mul(alpha, v)
                            );
  HF term_d3 = hf_mul(term_d3_inner, hf_mul(beta, z_omega_z));
  G1 d_3_s = g1_mul(&sigma_3_s, term_d3.value);

  // d_s = d_1_s + d_2_s - d_3_s
  G1 d_s = g1_add(&d_1_s, &d_2_s);
  G1 d_3_s_neg = g1_neg(&d_3_s);
  d_s = g1_add(&d_s, &d_3_s_neg);

  // ----------------------------------------------------
  // 10. Compute f_s
  // ----------------------------------------------------
  // We do scaled versions of t_lo_s, t_mid_s, t_hi_s
  G1 t_mid_s_scaled = g1_mul(&t_mid_s, hf_pow(z, n + 2).value);
  G1 t_hi_s_scaled = g1_mul(&t_hi_s, hf_pow(z, 2*n + 4).value);

  G1 t_lo_mid_s_scaled = g1_add(&t_lo_s, &t_mid_s_scaled);
  G1 t_hi_s_scaled_d_s = g1_add(&t_hi_s_scaled, &d_s);

  // a_s^(v^2), b_s^(v^3), c_s^(v^4), sigma1_s^(v^5), sigma2_s^(v^6)
  G1 a_s_pow_v2 = g1_mul(&a_s, hf_pow(v, 2).value);
  G1 b_s_pow_v3 = g1_mul(&b_s, hf_pow(v, 3).value);
  G1 c_s_pow_v4 = g1_mul(&c_s, hf_pow(v, 4).value);
  G1 sigma_1_s_pow_v5 = g1_mul(&sigma_1_s, hf_pow(v, 5).value);
  G1 sigma_2_s_pow_v6 = g1_mul(&sigma_2_s, hf_pow(v, 6).value);

  // combine them
  G1 temp_sum = g1_add(&a_s_pow_v2, &b_s_pow_v3);
  temp_sum = g1_add(&temp_sum, &c_s_pow_v4);
  temp_sum = g1_add(&temp_sum, &sigma_1_s_pow_v5);
  temp_sum = g1_add(&temp_sum, &sigma_2_s_pow_v6);

  // f_s = (t_lo_s + t_mid_s*z^(n+2)) + (t_hi_s*z^(2n+4) + d_s) + temp_sum
  G1 f_s = g1_add(&t_lo_mid_s_scaled, &t_hi_s_scaled_d_s);
  f_s = g1_add(&f_s, &temp_sum);

  // ----------------------------------------------------
  // 11. Compute e_s
  // ----------------------------------------------------
  // Build a polynomial "1" to do srs_eval_at_s
  HF one_val = hf_new(1);
  POLY one_poly = poly_new(&one_val, 1);
  G1 g1_s = srs_eval_at_s(&plonk->srs, &one_poly);
  poly_free(&one_poly);

  // e_scalar = t_z + v*r_z + v^2*a_z + v^3*b_z + v^4*c_z + v^5*s_sigma_1_z + v^6*s_sigma_2_z + u*z_omega_z
  HF e_scalar = hf_add(
      t_z,
      hf_add(
          hf_mul(v, r_z),
          hf_add(
              hf_mul(hf_pow(v, 2), a_z),
              hf_add(
                  hf_mul(hf_pow(v, 3), b_z),
                  hf_add(
                      hf_mul(hf_pow(v, 4), c_z),
                      hf_add(
                          hf_mul(hf_pow(v, 5), s_sigma_1_z),
                          hf_add(
                              hf_mul(hf_pow(v, 6), s_sigma_2_z),
                              hf_mul(u, z_omega_z)
                                 )
                             )
                         )
                     )
                 )
             )
                       );
  G1 e_s = g1_mul(&g1_s, e_scalar.value);

  // ----------------------------------------------------
  // 12. Batch validate pairings
  // ----------------------------------------------------
  // e_1 = pairing( (w_z_s + w_z_omega_s*u),  g2_s )
  // e_2 = pairing( (z*w_z_s + z*u*omega*w_z_omega_s) + (f_s - e_s), g2_1 )
  G1 w_z_omega_s_u = g1_mul(&w_z_omega_s, u.value);
  G1 e_1_q1 = g1_add(&w_z_s, &w_z_omega_s_u);
  G2 e_1_q2 = plonk->srs.g2_s;

  G1 temp1 = g1_mul(&w_z_s, z.value);
  G1 temp2 = g1_mul(&w_z_omega_s, hf_mul(u, hf_mul(z, omega)).value);
  G1 temp3 = g1_add(&temp1, &temp2);

  G1 e_s_neg = g1_neg(&e_s);
  G1 temp4 = g1_add(&f_s, &e_s_neg);
  G1 e_2_q1 = g1_add(&temp3, &temp4);
  G2 e_2_q2 = plonk->srs.g2_1;

  // Now do the pairing checks
  GTP e_1 = pairing(&e_1_q1, &e_1_q2);
  GTP e_2 = pairing(&e_2_q1, &e_2_q2);

  // If e_1 == e_2, success = true; else false
  success = gtp_equal(&e_1, &e_2);

cleanup:
  // ----------------------------------------------------
  // Free any resources still allocated
  // ----------------------------------------------------
  // sigma arrays (safe to free even if NULL)
  if (sigma_1) { free(sigma_1); sigma_1 = NULL; }
  if (sigma_2) { free(sigma_2); sigma_2 = NULL; }
  if (sigma_3) { free(sigma_3); sigma_3 = NULL; }

  // Polynomials
  poly_free(&q_m_x);
  poly_free(&q_l_x);
  poly_free(&q_r_x);
  poly_free(&q_o_x);
  poly_free(&q_c_x);
  poly_free(&s_sigma_1_x);
  poly_free(&s_sigma_2_x);
  poly_free(&s_sigma_3_x);

  return success;
}

#endif // PLONK_H
