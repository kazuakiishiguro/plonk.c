#ifndef PLONK_H
#define PLONK_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "constraints.h"
#include "matrix.h"
#include "pairing.h"
#include "poly.h"
#include "srs.h"
#include "hf.h"

#define OMEGA_VALUE 3 // example value; should be a primitive roor of unity in the field
#define K1_VALUE    2 // should not be in H
#define K2_VALUE    4 // should not be in H or K1 * H

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
  HF z_omega_z ;
} PROOF;

typedef struct {
  SRS srs;
  HF *h;          // array of hf elements
  size_t h_len;
  HF *k1_h;       // array of k1 * h elements
  HF *k2_h;       // array of k2 * h elements
  POLY z_h_x;     // polynomial z_h_x
} PLONK;

PLONK plonk_new(SRS srs, size_t n) {
  PLONK plonk;
  plonk.srs = srs;
  plonk.h_len = n + 1;

  // init
  HF OMEGA = hf_new(OMEGA_VALUE);
  HF K1 = hf_new(K1_VALUE);
  HF K2 = hf_new(K2_VALUE);

  // compute h = [omega^0, omega^1, ..., omega^n]
  plonk.h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < plonk.h_len; i++) {
       plonk.h[i] = hf_pow(OMEGA, i);
  }

  // check that K1 and K2 are not in H
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.h[i], K1) || hf_equal(plonk.h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }

  // compute h1_h: h1_h[i] = K1 * h[i]
  plonk.k1_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.k1_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k1_h[i] = hf_mul(K1, plonk.h[i]);
  }

  // compute k2_h: k2_h[i] = K2 * h[i]
  plonk.k2_h = (HF *)malloc(plonk.h_len * sizeof(HF));
  if (!plonk.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k2_h[i] = hf_mul(K2, plonk.h[i]);
  }

  // build h_pows MATRIX (vandermonde matrix) and compute its inverse
  MATRIX h_pows = matrix_zero(plonk.h_len, plonk.h_len);
  for (size_t r = 0; r < h_pows.m; r++) {
    HF h_pow = hf_one();
    for (size_t c = 0; c < h_pows.n; c++) {
      matrix_set(&h_pows, r, c, h_pow);
      h_pow = hf_mul(h_pow, plonk.h[r]);
    }
  }

  matrix_free(&h_pows);

  // compute z_h_x = Poly::z(h)
  plonk.z_h_x = poly_z(plonk.h, plonk.h_len);;

  return plonk;
}

void plonk_free(PLONK *plonk) {
  srs_free(&plonk->srs);
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
  poly_free(&plonk->z_h_x);
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

  // build vandermonde matrix over GF
  MATRIX h_pows = matrix_zero(plonk->h_len, plonk->h_len);
  for (size_t i = 0; i < plonk->h_len; i++) {
    HF x = plonk->h[i];  // roots of unity in HF
    HF x_pow = hf_one();
    for (size_t j = 0; j < plonk->h_len; j++) {
      matrix_set(&h_pows, i, j, x_pow); // set h_pows[i][j] = x_pow
      x_pow = hf_mul(x_pow, x); // x_pow *= x
    }
  }

  // invert the vanderomnde matrix over GF
  MATRIX h_pows_inv = matrix_inv(&h_pows); // ensure matrix_inverse works over GF
  matrix_free(&h_pows);

  // convert values to a column matrix over HF
  MATRIX vec = matrix_zero(len, 1);
  for (size_t i = 0; i < len; i++) {
    matrix_set(&vec, i , 0, values[i]);
  }

  // multiply h_pows_inv * vec over GF
  MATRIX res = matrix_mul(&h_pows_inv, &vec);
  matrix_free(&h_pows_inv);
  matrix_free(&vec);

  // convert result matrix to polynomial coefficients over HF
  HF *coeffs = (HF *)malloc(res.m * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }

  POLY result = poly_new(coeffs, res.m);
  free(coeffs);
  matrix_free(&res);

  return result;
}

PROOF plonk_prove(
    PLONK *plonk,
    CONSTRAINTS *constraints,
    ASSIGNMENTS *assignments,
    CHALLENGE *ch,
    HF rand[9]
                  ) {
  // step 1: Check that the constraints satisfy the assignments
  assert(constraints_satisfy(constraints, assignments));

  // extract challenges
  HF alpha = ch->alpha;
  HF beta = ch->beta;
  HF gamma = ch->gamma;
  HF z = ch->z;
  HF v = ch->v;

  // constants
  HF omega = hf_new(OMEGA_VALUE);
  HF k1 = hf_new(K1_VALUE);
  HF k2 = hf_new(K2_VALUE);

  size_t n = constraints->num_constraints;

  // step 2: create sigmas
  HF *sigma_1 = (HF *)malloc(n * sizeof(HF));
  HF *sigma_2 = (HF *)malloc(n * sizeof(HF));
  HF *sigma_3 = (HF *)malloc(n * sizeof(HF));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_prove\n");
    exit(EXIT_FAILURE);
  }

  copy_constraints_to_roots(plonk, constraints->c_a, n, sigma_1);
  copy_constraints_to_roots(plonk, constraints->c_b, n, sigma_2);
  copy_constraints_to_roots(plonk, constraints->c_c, n, sigma_3);

  // step 3: create polynomials

  // Pad a->a to length plonk->h_len
  HF *a_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!a_padded) {
    fprintf(stderr, "Memory allocation failed for a_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    a_padded[i] = f17(assignments->a[i].value);
  }

  POLY f_a_x = interpolate_at_h(plonk, a_padded, plonk->h_len);
  free(a_padded);

  HF *b_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!b_padded) {
    fprintf(stderr, "Memory allocation failed for b_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    b_padded[i] = f17(assignments->b[i].value);
  }

  POLY f_b_x = interpolate_at_h(plonk, b_padded, plonk->h_len);
  free(b_padded);

  HF *c_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!c_padded) {
    fprintf(stderr, "Memory allocation failed for c_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    c_padded[i] = f17(assignments->c[i].value);
  }
  POLY f_c_x = interpolate_at_h(plonk, c_padded, plonk->h_len);
  free(c_padded);

  HF *q_o_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *q_m_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *q_l_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *q_r_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *q_c_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!q_o_padded || !q_m_padded || !q_l_padded || !q_r_padded || !q_c_padded) {
       fprintf(stderr, "Memory allocation failed for padded constraint arrays\n");
       exit(EXIT_FAILURE);
  }
  memcpy(q_o_padded, constraints->q_o, n * sizeof(HF));
  memcpy(q_m_padded, constraints->q_m, n * sizeof(HF));
  memcpy(q_l_padded, constraints->q_l, n * sizeof(HF));
  memcpy(q_r_padded, constraints->q_r, n * sizeof(HF));
  memcpy(q_c_padded, constraints->q_c, n * sizeof(HF));

  POLY q_o_x = interpolate_at_h(plonk, q_o_padded, plonk->h_len);
  POLY q_m_x = interpolate_at_h(plonk, q_m_padded, plonk->h_len);
  POLY q_l_x = interpolate_at_h(plonk, q_l_padded, plonk->h_len);
  POLY q_r_x = interpolate_at_h(plonk, q_r_padded, plonk->h_len);
  POLY q_c_x = interpolate_at_h(plonk, q_c_padded, plonk->h_len);

  free(q_o_padded);
  free(q_m_padded);
  free(q_l_padded);
  free(q_r_padded);
  free(q_c_padded);

  // Similarly, pad sigma_1, sigma_2, sigma_3 before interpolating
  HF *sigma_1_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *sigma_2_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  HF *sigma_3_padded = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!sigma_1_padded || !sigma_2_padded || !sigma_3_padded) {
    fprintf(stderr, "Memory allocation failed for sigma_padded arrays\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    sigma_1_padded[i] = sigma_1[i];
    sigma_2_padded[i] = sigma_2[i];
    sigma_3_padded[i] = sigma_3[i];
  }

  POLY s_sigma_1 = interpolate_at_h(plonk, sigma_1_padded, plonk->h_len);
  POLY s_sigma_2 = interpolate_at_h(plonk, sigma_2_padded, plonk->h_len);
  POLY s_sigma_3 = interpolate_at_h(plonk, sigma_3_padded, plonk->h_len);
  free(sigma_1_padded);
  free(sigma_2_padded);
  free(sigma_3_padded);

  // step 4: round 1 - evaluate a(x), b(x), c(x) at s
  HF b1 = rand[0], b2 = rand[1], b3 = rand[2], b4 = rand[3], b5 = rand[4], b6 = rand[5];

  // a_x = (b2 + b1*x) * z_h_x + f_a_x
  HF a_blinding_coeffs[] = {b2, b1};
  POLY a_blinding_poly  = poly_new(a_blinding_coeffs, 2);
  POLY a_x_blinded = poly_mul(&a_blinding_poly, &plonk->z_h_x);
  POLY a_x = poly_add(&a_x_blinded, &f_a_x);

  // b_x = (b4 + b3*x) * z_h_x + f_b_x
  HF b_blinding_coeffs[] = {b4, b3};
  POLY b_blinding_poly  = poly_new(b_blinding_coeffs, 2);
  POLY b_x_blinded = poly_mul(&b_blinding_poly, &plonk->z_h_x);
  POLY b_x = poly_add(&b_x_blinded, &f_b_x);

  // c_x = (b6 + b5*x) * z_h_x + f_c_x
  HF c_blinding_coeffs[] = {b6, b5};
  POLY c_blinding_poly  = poly_new(c_blinding_coeffs, 2);
  POLY c_x_blinded = poly_mul(&c_blinding_poly, &plonk->z_h_x);
  POLY c_x = poly_add(&c_x_blinded, &f_c_x);

  // evaluate at s
  G1 a_s = srs_eval_at_s(&plonk->srs, &a_x);
  G1 b_s = srs_eval_at_s(&plonk->srs, &b_x);
  G1 c_s = srs_eval_at_s(&plonk->srs, &c_x);

  // cleaning up
  poly_free(&a_blinding_poly);
  poly_free(&a_x_blinded);
  poly_free(&b_blinding_poly);
  poly_free(&b_x_blinded);
  poly_free(&c_blinding_poly);
  poly_free(&c_x_blinded);

  // step 5: round 2 - evaluate the accumultor vector polynomial at s
  HF b7 = rand[6], b8 = rand[7], b9 = rand[8];

  // initialize accumulator vector with 1
  HF *acc = (HF *)malloc(plonk->h_len *sizeof(HF));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    exit(EXIT_FAILURE);
  }
  acc[0] = hf_new(1);

  for (size_t i = 1; i <= n; i++) {
    // assignments at positoin i - 1
    HF as_a = assignments->a[i - 1];
    HF as_b = assignments->b[i - 1];
    HF as_c = assignments->c[i - 1];

    // omega_pow = omega^(i-1)
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
    HF s_sigma_1_eval = poly_eval(&s_sigma_1, f17(omega_pow.value));
    HF s_sigma_2_eval = poly_eval(&s_sigma_2, f17(omega_pow.value));
    HF s_sigma_3_eval = poly_eval(&s_sigma_3, f17(omega_pow.value));

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

  // pad the rest with zeros
  for (size_t i = n + 1; i < plonk->h_len; i++)
    acc[i] = hf_new(0);

  // interpolate accumualtor vector
  POLY acc_x = interpolate_at_h(plonk, acc, plonk->h_len);
  free(acc);

  // check that acc_x evaluated at omega^n is 1
  HF omega_n = hf_pow(omega, n);
  printf("omega_n = %u\n", omega_n.value);
  HF acc_x_eval = poly_eval(&acc_x, omega_n);
  printf("acc_x_eval = %u\n", acc_x_eval.value);
  printf("Evaluating acc_x at x = %u\n", omega_n.value);
  printf("acc_x coefficients:\n");
  for (size_t i = 0; i < acc_x.len; i++) {
    printf("Coefficient of x^%zu: %u\n", i, acc_x.coeffs[i].value);
  }
  printf("acc values:\n");
  for (size_t i = 0; i < plonk->h_len; i++) {
    printf("acc[%zu] = %u\n", i, acc[i].value);
  }

  assert(hf_equal(acc_x_eval, hf_one()));

  // create z_x polynomial
  HF z_blinding_coeffs[] = {b9, b8, b7};
  POLY z_blinding_poly = poly_new(z_blinding_coeffs, 3);
  POLY z_blinded = poly_mul(&z_blinding_poly, &plonk->z_h_x);
  POLY z_x = poly_add(&z_blinded, &acc_x);

  // evaluate z_x at s
  G1 z_s = srs_eval_at_s(&plonk->srs, &z_x);

  // clean up
  poly_free(&z_blinding_poly);
  poly_free(&z_blinded);

  // step 6: Compute the quotient polynomial t(x)
  // compute Lagrange polynomial L1(x)
  HF *lagrange_vector = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed for lagrange_vector\n");
    exit(EXIT_FAILURE);
  }
  lagrange_vector[0] = hf_one(); // L1(omega^0) = 1, rest are 0
  POLY l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);
  free(lagrange_vector);

  // compute p_i_x (public input polynomial)
  // assuming no public inputs for simplicity
  POLY p_i_x = poly_zero();

  // compute helper polynomials
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
  POLY tmp = poly_new(&gamma, 1);
  POLY a_x_beta = poly_add(&a_x, &tmp);
  tmp = poly_new(&omega, 1);
  tmp = poly_scale(&tmp, beta);
  a_x_beta = poly_add(&a_x_beta, &tmp);
  tmp = poly_new(&gamma, 1);
  POLY b_x_beta_k1 = poly_add(&b_x, &tmp);
  HF beta_k1 = hf_mul(beta, k1);
  POLY beta_k1_1 = poly_new(&beta_k1, 1);
  POLY beta_k1_1_1 = poly_scale(&beta_k1_1, hf_one());
  b_x_beta_k1 = poly_add(&b_x_beta_k1, &beta_k1_1_1);
  POLY c_x_beta_k2 = poly_add(&c_x, &tmp);
  HF beta_k2 = hf_mul(beta, k2);
  POLY beta_k2_1 = poly_new(&beta_k2, 1);
  POLY beta_k2_1_1 = poly_scale(&beta_k2_1, hf_one());
  c_x_beta_k2 = poly_add(&c_x_beta_k2, &beta_k2_1_1);
  POLY temp = poly_mul(&a_x_beta, &b_x_beta_k1);
  POLY t_2_z_h = poly_mul(&temp, &c_x_beta_k2);
  poly_free(&temp);
  poly_free(&a_x_beta);
  poly_free(&beta_k1_1);
  poly_free(&beta_k1_1_1);
  poly_free(&b_x_beta_k1);
  poly_free(&beta_k2_1);
  poly_free(&beta_k2_1_1);
  poly_free(&c_x_beta_k2);
  t_2_z_h = poly_scale(&t_2_z_h, alpha);
  t_2_z_h = poly_mul(&t_2_z_h, &z_x);

  // compute t_3_z_h
  POLY s_sigma_1_beta = poly_scale(&s_sigma_1, beta);
  tmp = poly_new(&gamma, 1);
  s_sigma_1_beta = poly_add(&s_sigma_1_beta, &tmp);
  POLY s_sigma_2_beta = poly_scale(&s_sigma_2, beta);
  s_sigma_2_beta = poly_add(&s_sigma_2_beta, &tmp);
  POLY s_sigma_3_beta = poly_scale(&s_sigma_3, beta);
  s_sigma_3_beta = poly_add(&s_sigma_3_beta, &tmp);
  poly_free(&tmp);

  POLY temp1 = poly_add(&a_x, &s_sigma_1_beta);
  POLY temp2 = poly_add(&b_x, &s_sigma_2_beta);
  POLY temp3 = poly_add(&c_x, &s_sigma_3_beta);
  poly_free(&s_sigma_1_beta);
  poly_free(&s_sigma_2_beta);
  poly_free(&s_sigma_3_beta);

  POLY temp4 = poly_mul(&temp1, &temp2);
  POLY temp5 = poly_mul(&temp4, &temp3);
  poly_free(&temp1);
  poly_free(&temp2);
  poly_free(&temp3);
  poly_free(&temp4);

  POLY z_omega_x = poly_shift(&z_x, 1); // Shift z_x by multiplying x^1
  POLY t_3_z_h = poly_mul(&temp5, &z_omega_x);
  poly_free(&temp5);
  poly_free(&z_omega_x);
  t_3_z_h = poly_scale(&t_3_z_h, alpha);

  // compute t_4_z_h
  HF one = hf_one();
  POLY one_1 = poly_new(&one, 1);
  POLY one_1_neg = poly_negate(&one_1);
  POLY z_x_minus_one = poly_add(&z_x, &one_1_neg);
  POLY t_4_z_h = poly_scale(&z_x_minus_one, hf_pow(alpha, 2));
  poly_free(&one_1);
  poly_free(&one_1_neg);
  poly_free(&z_x_minus_one);
  t_4_z_h = poly_mul(&t_4_z_h, &l_1_x);

  // compute t(x)
  POLY t_x_numerator = poly_add(&t_1_z_h, &t_2_z_h);
  poly_free(&t_1_z_h);
  poly_free(&t_2_z_h);
  t_x_numerator = poly_sub(&t_x_numerator, &t_3_z_h);
  poly_free(&t_3_z_h);
  t_x_numerator = poly_add(&t_x_numerator, &t_4_z_h);
  poly_free(&t_4_z_h);

  // divide t_x_numerator by z_h_x to get t_x
  POLY t_x;
  POLY remainder;
  poly_divide(&t_x_numerator, &plonk->z_h_x, &t_x, &remainder);
  poly_free(&t_x_numerator);
  if (!poly_is_zero(&remainder)) {
    fprintf(stderr, "Non-zero remainder in t(x) division\n");
    exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // evaluate t_x at z to get t_z
  HF t_z = poly_eval(&t_x, z);

  // step 7: Split t(x) into t_lo_x, t_mid_x, t_hi_x
  size_t degree = t_x.len;
  size_t n_plus_2 = n + 2;
  size_t part_size = n_plus_2;
  size_t total_parts = 3;

  // ensure t_x has enough coefficients
  if (degree < 3 * part_size) {
    // pad t_x with zeros if necessary
    HF *padded_coeffs = (HF *)calloc(3 * part_size, sizeof(HF));
    if (!padded_coeffs) {
      fprintf(stderr, "Memory allocation failed in padding t_x\n");
      exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < degree; i++) {
      padded_coeffs[i] = t_x.coeffs[i];
    }
    poly_free(&t_x);
    t_x = poly_new(padded_coeffs, 3 * part_size);
    free(padded_coeffs);
  }

  POLY t_lo_x = poly_slice(&t_x, 0, part_size);
  POLY t_mid_x = poly_slice(&t_x, part_size, 2 * part_size);
  POLY t_hi_x = poly_slice(&t_x, 2 * part_size, degree);
  poly_free(&t_x);

  // step 8: Evaluate t_lo_x, t_mid_x, t_hi_x at s
  G1 t_lo_s = srs_eval_at_s(&plonk->srs, &t_lo_x);
  G1 t_mid_s = srs_eval_at_s(&plonk->srs, &t_mid_x);
  G1 t_hi_s = srs_eval_at_s(&plonk->srs, &t_hi_x);

  // clean up t polynomials
  poly_free(&t_lo_x);
  poly_free(&t_mid_x);
  poly_free(&t_hi_x);

  // Step 9: Compute r(x) and evaluate it at z
  // Evaluate polynomials at z
  HF a_z = poly_eval(&a_x, z);
  HF b_z = poly_eval(&b_x, z);
  HF c_z = poly_eval(&c_x, z);
  HF s_sigma_1_z = poly_eval(&s_sigma_1, z);
  HF s_sigma_2_z = poly_eval(&s_sigma_2, z);
  HF z_omega_z = poly_eval(&z_x, hf_mul(z, omega));

  // Compute r_1_x
  POLY q_m_x_a_z_b_z = poly_scale(&q_m_x, hf_mul(a_z, b_z));
  POLY q_l_x_a_z = poly_scale(&q_l_x, a_z);
  POLY q_r_x_b_z = poly_scale(&q_r_x, b_z);
  POLY q_o_x_c_z = poly_scale(&q_o_x, c_z);
  POLY sum_r1 = poly_add(&q_m_x_a_z_b_z, &q_l_x_a_z);
  poly_free(&q_m_x_a_z_b_z);
  poly_free(&q_l_x_a_z);
  sum_r1 = poly_add(&sum_r1, &q_r_x_b_z);
  poly_free(&q_r_x_b_z);
  sum_r1 = poly_add(&sum_r1, &q_o_x_c_z);
  poly_free(&q_o_x_c_z);
  POLY r_1_x = poly_add(&sum_r1, &q_c_x);
  poly_free(&sum_r1);

  // Compute r_2_x
  HF beta_z = hf_mul(beta, z);
  HF gamma_beta_z = hf_add(gamma, beta_z);
  HF beta_k1_z = hf_mul(beta, hf_mul(k1, z));
  HF gamma_beta_k1_z = hf_add(gamma, beta_k1_z);
  HF beta_k2_z = hf_mul(beta, hf_mul(k2, z));
  HF gamma_beta_k2_z = hf_add(gamma, beta_k2_z);

  HF term1 = hf_add(a_z, gamma_beta_z);
  HF term2 = hf_add(b_z, gamma_beta_k1_z);
  HF term3 = hf_add(c_z, gamma_beta_k2_z);

  HF numerator = hf_mul(hf_mul(term1, term2), term3);
  POLY r_2_x = poly_scale(&z_x, hf_mul(numerator, alpha));

  // Compute r_3_x
  HF beta_s_sigma_1_z = hf_mul(beta, s_sigma_1_z);
  HF gamma_beta_s_sigma_1_z = hf_add(gamma, beta_s_sigma_1_z);
  HF beta_s_sigma_2_z = hf_mul(beta, s_sigma_2_z);
  HF gamma_beta_s_sigma_2_z = hf_add(gamma, beta_s_sigma_2_z);

  HF term4 = hf_add(a_z, gamma_beta_s_sigma_1_z);
  HF term5 = hf_add(b_z, gamma_beta_s_sigma_2_z);
  HF temp_alpha = hf_mul(alpha, alpha);

  HF s_sigma_3_z = poly_eval(&s_sigma_3, z);
  HF beta_s_sigma_3_z = hf_mul(beta, s_sigma_3_z);
  HF gamma_beta_s_sigma_3_z = hf_add(gamma, beta_s_sigma_3_z);

  HF term6 = hf_mul(z_omega_z, beta_s_sigma_3_z);

  HF numerator_r3 = hf_mul(hf_mul(term4, term5), term6);
  POLY r_3_x = poly_scale(&z_x, hf_mul(numerator_r3, alpha));

  // Compute r_4_x
  HF l_1_z = poly_eval(&l_1_x, z);
  POLY r_4_x = poly_scale(&z_x, hf_mul(hf_mul(temp_alpha, l_1_z), hf_new(1)));

  // Compute r_x
  POLY sum_r = poly_add(&r_1_x, &r_2_x);
  poly_free(&r_1_x);
  poly_free(&r_2_x);
  sum_r = poly_add(&sum_r, &r_3_x);
  poly_free(&r_3_x);
  sum_r = poly_add(&sum_r, &r_4_x);
  poly_free(&r_4_x);

  POLY r_x = sum_r;

  // Evaluate r_x at z to get r_z
  HF r_z = poly_eval(&r_x, z);

  // Step 10: Compute opening prf w_z_x and w_z_omega_x
  // Compute t_lo_x + t_mid_x * z^{n+2} + t_hi_x * z^{2n+4}
  HF z_pow_n_plus_2 = hf_pow(z, n + 2);
  HF z_pow_2n_plus_4 = hf_pow(z, 2 * n + 4);

  POLY t_mid_x_z = poly_scale(&t_mid_x, z_pow_n_plus_2);
  POLY t_hi_x_z = poly_scale(&t_hi_x, z_pow_2n_plus_4);

  POLY t_sum = poly_add(&t_lo_x, &t_mid_x_z);
  poly_free(&t_mid_x_z);
  t_sum = poly_add(&t_sum, &t_hi_x_z);
  poly_free(&t_hi_x_z);

  // Compute w_z_x numerator
  POLY t_z_poly = poly_new(&t_z, 1); // Polynomial representing the constant t_z
  POLY w_z_x_num = poly_sub(&t_sum, &t_z_poly);
  poly_free(&t_sum);
  poly_free(&t_z_poly);

  // Compute (r_x - r_z) * v
  POLY r_z_1 = poly_new(&r_z, 1);
  POLY r_x_minus_r_z = poly_sub(&r_x, &r_z_1);
  POLY r_x_minus_r_z_v = poly_scale(&r_x_minus_r_z, v);
  poly_free(&r_z_1);
  poly_free(&r_x_minus_r_z);

  // Compute (a_x - a_z) * v^2
  POLY a_z_1 = poly_new(&a_z, 1);
  POLY a_x_minus_a_z = poly_sub(&a_x, &a_z_1);
  HF v_squared = hf_mul(v, v);
  POLY a_x_minus_a_z_v2 = poly_scale(&a_x_minus_a_z, v_squared);
  poly_free(&a_z_1);
  poly_free(&a_x_minus_a_z);

  // Similarly compute (b_x - b_z) * v^3 and (c_x - c_z) * v^4
  HF v_cubed = hf_mul(v_squared, v);
  POLY b_z_1 = poly_new(&b_z, 1);
  POLY b_x_minus_b_z = poly_sub(&b_x, &b_z_1);
  POLY b_x_minus_b_z_v3 = poly_scale(&b_x_minus_b_z, v_cubed);
  poly_free(&b_z_1);
  poly_free(&b_x_minus_b_z);

  HF v_fourth = hf_mul(v_cubed, v);
  POLY c_z_1 = poly_new(&c_z, 1);
  POLY c_x_minus_c_z = poly_sub(&c_x, &c_z_1);
  POLY c_x_minus_c_z_v4 = poly_scale(&c_x_minus_c_z, v_fourth);
  poly_free(&c_z_1);
  poly_free(&c_x_minus_c_z);

  // Similarly compute terms for s_sigma_1 and s_sigma_2
  HF v_fifth = hf_mul(v_fourth, v);
  POLY s_sigma_1_z_1 = poly_new(&s_sigma_1_z, 1);
  POLY s_sigma_1_x_minus_s_sigma_1_z = poly_sub(&s_sigma_1, &s_sigma_1_z_1);
  POLY s_sigma_1_x_minus_s_sigma_1_z_v5 = poly_scale(&s_sigma_1_x_minus_s_sigma_1_z, v_fifth);
  poly_free(&s_sigma_1_z_1);
  poly_free(&s_sigma_1_x_minus_s_sigma_1_z);

  HF v_sixth = hf_mul(v_fifth, v);
  POLY s_sigma_2_z_1 = poly_new(&s_sigma_2_z, 1);
  POLY s_sigma_2_x_minus_s_sigma_2_z = poly_sub(&s_sigma_2, &s_sigma_2_z_1);
  POLY s_sigma_2_x_minus_s_sigma_2_z_v6 = poly_scale(&s_sigma_2_x_minus_s_sigma_2_z, v_sixth);
  poly_free(&s_sigma_2_z_1);
  poly_free(&s_sigma_2_x_minus_s_sigma_2_z);

  // Sum all terms to compute w_z_x numerator
  POLY w_z_x_num_total = poly_add(&w_z_x_num, &r_x_minus_r_z_v);
  poly_free(&w_z_x_num);
  poly_free(&r_x_minus_r_z_v);
  w_z_x_num_total = poly_add(&w_z_x_num_total, &a_x_minus_a_z_v2);
  poly_free(&a_x_minus_a_z_v2);
  w_z_x_num_total = poly_add(&w_z_x_num_total, &b_x_minus_b_z_v3);
  poly_free(&b_x_minus_b_z_v3);
  w_z_x_num_total = poly_add(&w_z_x_num_total, &c_x_minus_c_z_v4);
  poly_free(&c_x_minus_c_z_v4);
  w_z_x_num_total = poly_add(&w_z_x_num_total, &s_sigma_1_x_minus_s_sigma_1_z_v5);
  poly_free(&s_sigma_1_x_minus_s_sigma_1_z_v5);
  w_z_x_num_total = poly_add(&w_z_x_num_total, &s_sigma_2_x_minus_s_sigma_2_z_v6);
  poly_free(&s_sigma_2_x_minus_s_sigma_2_z_v6);

  // Divide w_z_x numerator by (x - z)
  POLY x_minus_z = poly_new((HF[]){hf_neg(z), hf_new(1)}, 2);
  POLY w_z_x;
  POLY rem;
  poly_divide(&w_z_x_num_total, &x_minus_z, &w_z_x, &rem);
  poly_free(&w_z_x_num_total);
  poly_free(&x_minus_z);
  if (!poly_is_zero(&remainder)) {
       fprintf(stderr, "Non-zero remainder in w_z_x division\n");
       exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // Compute w_z_s = srs_eval_at_s(w_z_x)
  G1 w_z_s = srs_eval_at_s(&plonk->srs, &w_z_x);

  // Compute w_z_omega_x
  POLY z_omega_z_1 = poly_new(&z_omega_z, 1);
  POLY z_x_minus_z_omega_z = poly_sub(&z_x, &z_omega_z_1);
  HF z_omega = hf_mul(z, omega);
  POLY x_minus_z_omega = poly_new((HF[]){hf_neg(z_omega), hf_one()}, 2);
  POLY w_z_omega_x;
  poly_divide(&z_x_minus_z_omega_z, &x_minus_z_omega, &w_z_omega_x, &remainder);
  poly_free(&z_omega_z_1);
  poly_free(&z_x_minus_z_omega_z);
  poly_free(&x_minus_z_omega);
  if (!poly_is_zero(&remainder)) {
       fprintf(stderr, "Non-zero remainder in w_z_omega_x division\n");
       exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // Compute w_z_omega_s = srs_eval_at_s(w_z_omega_x)
  G1 w_z_omega_s = srs_eval_at_s(&plonk->srs, &w_z_omega_x);

  // Clean up polynomials
  poly_free(&w_z_x);
  poly_free(&w_z_omega_x);

  // Step 11: Assemble the Proof
  PROOF prf = {
       .a_s = a_s,
       .b_s = b_s,
       .c_s = c_s,
       .z_s = z_s,
       .t_lo_s = t_lo_s,
       .t_mid_s = t_mid_s,
       .t_hi_s = t_hi_s,
       .w_z_s = w_z_s,
       .w_z_omega_s = w_z_omega_s,
       .a_z = a_z,
       .b_z = b_z,
       .c_z = c_z,
       .s_sigma_1_z = s_sigma_1_z,
       .s_sigma_2_z = s_sigma_2_z,
       .r_z = r_z,
       .z_omega_z = z_omega_z
  };

  // Clean up remaining polynomials
  poly_free(&a_x);
  poly_free(&b_x);
  poly_free(&c_x);
  poly_free(&z_x);
  poly_free(&r_x);
  poly_free(&l_1_x);
  poly_free(&s_sigma_1);
  poly_free(&s_sigma_2);
  poly_free(&s_sigma_3);
  poly_free(&q_l_x);
  poly_free(&q_r_x);
  poly_free(&q_o_x);
  poly_free(&q_m_x);
  poly_free(&q_c_x);

  return prf;
}

bool plonk_verify(
    PLONK *plonk,
    CONSTRAINTS *cons,
    PROOF *prf,
    CHALLENGE *ch,
    HF rand[1]
                  ) {
  // unpack prf
  G1 a_s = prf->a_s;
  G1 b_s = prf->b_s;
  G1 c_s = prf->c_s;
  G1 z_s = prf->z_s;
  G1 t_lo_s = prf->t_lo_s;
  G1 t_mid_s = prf->t_mid_s;
  G1 t_hi_s = prf->t_hi_s;
  G1 w_z_s = prf->w_z_s;
  G1 w_z_omega_s = prf->w_z_omega_s;
  HF a_z = prf->a_z;
  HF b_z = prf->b_z;
  HF c_z = prf->c_z;
  HF s_sigma_1_z = prf->s_sigma_1_z;
  HF s_sigma_2_z = prf->s_sigma_2_z;
  HF r_z = prf->r_z;
  HF z_omega_z = prf->z_omega_z;

  // unpack challenge
  HF alpha = ch->alpha;
  HF beta = ch->beta;
  HF gamma = ch->gamma;
  HF z = ch->z;
  HF v = ch->v;

  // constants
  HF omega = hf_new(OMEGA_VALUE);
  HF k1 = hf_new(K1_VALUE);
  HF k2 = hf_new(K2_VALUE);

  // step 1: Validate proof points in G1
  /* if (!g1_p_is_on_curve(&a_s) || */
  /*     !g1_p_is_on_curve(&b_s) || */
  /*     !g1_p_is_on_curve(&c_s) || */
  /*     !g1_p_is_on_curve(&z_s) || */
  /*     !g1_p_is_on_curve(&t_lo_s) || */
  /*     !g1_p_is_on_curve(&t_mid_s) || */
  /*     !g1_p_is_on_curve(&t_hi_s) || */
  /*     !g1_p_is_on_curve(&w_z_s) || */
  /*     !g1_p_is_on_curve(&w_z_omega_s)) { */
  /*   return false; */
  /* } */

  // step 2: Validate proof fields in HF
  /* if (!hf_in_field(a_z) || */
  /*     !hf_in_field(b_z) || */
  /*     !hf_in_field(c_z) || */
  /*     !hf_in_field(s_sigma_1_z) || */
  /*     !hf_in_field(s_sigma_2_z) || */
  /*     !hf_in_field(r_z) || */
  /*     !hf_in_field(z_omega_z)) { */
  /*   return false; */
  /* } */

  // step 3: No public inputs, nothing to do

  // step 4: Evaluate z_h at z
  HF z_h_z = poly_eval(&plonk->z_h_x, z);

  // step 5: Evaluate lagrange polynomial L1 at z
  HF *lagrange_vector = (HF *)calloc(plonk->h_len, sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    exit(EXIT_FAILURE);
  }

  lagrange_vector[0] = hf_new(1); // L1(omega^0) = 1
  POLY l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);
  free(lagrange_vector);
  HF l_1_z = poly_eval(&l_1_x, z);
  poly_free(&l_1_x);

  // step 6: No public inputs, nothing to do
  HF p_i_z = hf_new(0);

  // step 7: Compute quotient polynomial evaluation
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_1_z), gamma), a_z);
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_2_z), gamma), b_z);
  HF c_z_gamma = hf_add(c_z, gamma);
  HF l_1_z_alpha_2 = hf_mul(l_1_z, hf_pow(alpha, 2));

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

  HF t_z = hf_div(numerator, z_h_z);

  // step 8: Compute the first part of batched polynomial commitment
  // compute sigmas
  size_t n = cons->num_constraints;

  HF *sigma_1 = (HF *)malloc(n * sizeof(HF));
  HF *sigma_2 = (HF *)malloc(n * sizeof(HF));
  HF *sigma_3 = (HF *)malloc(n * sizeof(HF));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    exit(EXIT_FAILURE);
  }

  copy_constraints_to_roots(plonk, cons->c_a, n, sigma_1);
  copy_constraints_to_roots(plonk, cons->c_b, n, sigma_2);
  copy_constraints_to_roots(plonk, cons->c_c, n, sigma_3);

  // interpolate and evaluate polynomials at s
  POLY q_m_x = interpolate_at_h(plonk, cons->q_m, n);
  G1 q_m_s = srs_eval_at_s(&plonk->srs, &q_m_x);
  poly_free(&q_m_x);

  POLY q_l_x = interpolate_at_h(plonk, cons->q_l, n);
  G1 q_l_s = srs_eval_at_s(&plonk->srs, &q_l_x);
  poly_free(&q_l_x);

  POLY q_r_x = interpolate_at_h(plonk, cons->q_r, n);
  G1 q_r_s = srs_eval_at_s(&plonk->srs, &q_r_x);
  poly_free(&q_r_x);

  POLY q_o_x = interpolate_at_h(plonk, cons->q_o, n);
  G1 q_o_s = srs_eval_at_s(&plonk->srs, &q_o_x);
  poly_free(&q_o_x);

  POLY q_c_x = interpolate_at_h(plonk, cons->q_c, n);
  G1 q_c_s = srs_eval_at_s(&plonk->srs, &q_c_x);
  poly_free(&q_c_x);

  POLY s_sigma_1_x = interpolate_at_h(plonk, sigma_1, n);
  G1 sigma_1_s = srs_eval_at_s(&plonk->srs, &s_sigma_1_x);
  poly_free(&s_sigma_1_x);

  POLY s_sigma_2_x = interpolate_at_h(plonk, sigma_2, n);
  G1 sigma_2_s = srs_eval_at_s(&plonk->srs, &s_sigma_2_x);
  poly_free(&s_sigma_2_x);

  POLY s_sigma_3_x = interpolate_at_h(plonk, sigma_3, n);
  G1 sigma_3_s = srs_eval_at_s(&plonk->srs, &s_sigma_3_x);
  poly_free(&s_sigma_3_x);

  free(sigma_1);
  free(sigma_2);
  free(sigma_3);

  // Random scalar u
  HF u = rand[0];

  // Compute d_1_s
  G1 tmp_lhs = g1_mul(&q_r_s, hf_mul(b_z, v).value);
  G1 tmp_rhs = g1_mul(&q_o_s, hf_mul(c_z, v).value);
  G1 tmp_lhs2 = g1_mul(&q_l_s, hf_mul(a_z, v).value);
  G1 tmp_rhs2 = g1_add(&tmp_lhs, &tmp_rhs);
  G1 tmp_lhs3 = g1_mul(&q_m_s, hf_mul(hf_mul(a_z, b_z), v).value);
  G1 tmp_rhs3 = g1_add(&tmp_lhs2, &tmp_rhs2);
  G1 d_1_s = g1_add(&tmp_lhs3,&tmp_rhs3);
  G1 tmp_rhs4 = g1_mul(&q_c_s, v.value);
  d_1_s = g1_add(&d_1_s, &tmp_rhs4);

  // Compute d_2_s
  HF term_d2 = hf_add(
      hf_mul(
          hf_mul(
              hf_mul(
                  hf_add(a_z, hf_add(hf_mul(beta, z), gamma)),
                  hf_add(b_z, hf_add(hf_mul(beta, hf_mul(k1, z)), gamma))
                        ),
              hf_add(c_z, hf_add(hf_mul(beta, hf_mul(k2, z)), gamma))
                    ),
          hf_mul(alpha, v)
                ),
      hf_mul(l_1_z, hf_mul(hf_pow(alpha, 2), v))
                            );
  term_d2 = hf_add(term_d2, u);
  G1 d_2_s = g1_mul(&z_s, term_d2.value);

  // Compute d_3_s
  HF term_d3 = hf_mul(
      hf_mul(
          hf_mul(
              hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma)),
              hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma))
                    ),
          hf_mul(alpha, v)
                ),
      hf_mul(beta, z_omega_z)
                            );
  G1 d_3_s = g1_mul(&sigma_3_s, term_d3.value);

  // Compute d_s = d_1_s + d_2_s - d_3_s
  G1 d_s = g1_add(&d_1_s, &d_2_s);
  G1 d_3_s_neg = g1_neg(&d_3_s);
  d_s = g1_add(&d_s, &d_3_s_neg);

  // Step 9: Compute f_s
  G1 t_lo_s_scaled = t_lo_s;
  G1 t_mid_s_scaled = g1_mul(&t_mid_s, hf_pow(z, n + 2).value);
  G1 t_hi_s_scaled = g1_mul(&t_hi_s, hf_pow(z, 2 * n + 4).value);
  G1 t_lo_mid_s_scaled = g1_add(&t_lo_s_scaled, &t_mid_s_scaled);
  G1 t_hi_s_scaled_d_s = g1_add(&t_hi_s_scaled, &d_s);
  G1 a_s_pow_v_2 = g1_mul(&a_s, hf_pow(v, 2).value);
  G1 b_s_pow_v_3 = g1_mul(&b_s, hf_pow(v, 3).value);
  G1 c_s_pow_v_4 = g1_mul(&c_s, hf_pow(v, 4).value);
  G1 sigma_1_s_pow_v_5 = g1_mul(&sigma_1_s, hf_pow(v, 5).value);
  G1 sigma_2_s_pow_v_6 = g1_mul(&sigma_2_s, hf_pow(v, 6).value);
  G1 temp_sum = g1_add(&sigma_1_s_pow_v_5, &sigma_2_s_pow_v_6);
  temp_sum = g1_add(&c_s_pow_v_4, &temp_sum);
  temp_sum = g1_add(&b_s_pow_v_3, &temp_sum);
  temp_sum = g1_add(&a_s_pow_v_2, &temp_sum);

  G1 f_s = g1_add(&t_lo_mid_s_scaled, &t_hi_s_scaled_d_s);
  f_s = g1_add(&f_s, &temp_sum);

  // Step 10: Compute e_s
  POLY one_poly = poly_new(&(HF){ .value = 1 }, 1);
  G1 g1_s = srs_eval_at_s(&plonk->srs, &one_poly);
  poly_free(&one_poly);

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

  // step 11: Batch validate all equations
  G1 w_z_omega_s_u = g1_mul(&w_z_omega_s, u.value);
  G1 e_1_q1 = g1_add(&w_z_s, &w_z_omega_s_u);
  G2 e_1_q2 = plonk->srs.g2_s;

  G1 temp1 = g1_mul(&w_z_s, z.value);
  G1 temp2 = g1_mul(&w_z_omega_s, hf_mul(u, hf_mul(z, omega)).value);
  G1 temp3 = g1_add(&temp1, &temp2);

  // TODO: check it's validity for g1_p_sub
  G1 e_s_neg = g1_neg(&e_s);
  G1 temp4 = g1_add(&f_s, &e_s_neg);

  G1 e_2_q1 = g1_add(&temp3, &temp4);
  G2 e_2_q2 = plonk->srs.g2_1;

  // compute pairings
  GTP e_1 = pairing(&e_1_q1, &e_1_q2);
  GTP e_2 = pairing(&e_2_q1, &e_2_q2);

  // check if e_1 == e_2
  if (gtp_equal(&e_1, &e_2)) {
    return true;
  } else {
    return false;
  }
}

#endif // PLONK_H
