#ifndef PLONK_H
#define PLONK_H

#include <assert.h>
#include <stdlib.h>
#include "constraints.h"
#include "matrix.h"
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
  // step 1: check that the constraints satisfies the assignments
  assert(constraints_satisfy(constraints, assignments));

  // extract challenges
  HF alpha = challenge->alpha;
  HF beta = challenge->beta;
  HF gamma = challenge->gamma;
  HF z = challenge->z;
  HF v = challenge->v;

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

  // step 3: create a set of polynomials, those polynomials evaluates at roots of unity for all the components of the plonk circuit
  //
  //   (a,b,c)                  : assignments
  //   (o,m,l,r,c)              : gate constraints
  //   (sigma1, sigma2, sigma3) : copy constraints
  //
  // TODO: check if these elementes need to be padded
  POLY f_a_x = interpolate_at_h(plonk, assignments->a, plonk->h_len);
  POLY f_b_x = interpolate_at_h(plonk, assignments->b, plonk->h_len);
  POLY f_c_x = interpolate_at_h(plonk, assignments->c, plonk->h_len);
  POLY q_o_x = interpolate_at_h(plonk, constraints->q_o, plonk->h_len);
  POLY q_m_x = interpolate_at_h(plonk, constraints->q_m, plonk->h_len);
  POLY q_l_x = interpolate_at_h(plonk, constraints->q_l, plonk->h_len);
  POLY q_r_x = interpolate_at_h(plonk, constraints->q_r, plonk->h_len);
  POLY q_c_x = interpolate_at_h(plonk, constraints->q_c, plonk->h_len);
  POLY s_sigma_1 = interpolate_at_h(plonk, sigma_1, plonk->h_len);
  POLY s_sigma_2 = interpolate_at_h(plonk, sigma_2, plonk->h_len);
  POLY s_sigma_3 = interpolate_at_h(plonk, sigma_3, plonk->h_len);

  // step 4: round 1 - eval a(x), b(x), c(x) at s
  /* ------------------------------------------------------- */
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

  // output of first step
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

  // round 2 - eval the accumulator vector polynomial at s
  /* ------------------------------------------------------- */
  // check https://vitalik.ca/general/2019/09/22/plonk.html

  // initialize the accumulator vector with 1
  //
  //   beta   blinding against prover addition manipulation
  //   gamma  blinding against prover multiplication manipulation
  //   k1,k2  coset independance between a,b,c
  HF *acc = (HF *)malloc(n * sizeof(HF));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    exit(EXIT_FAILURE);
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
    HF s_sigma_1_eval = poly_eval(&s_sigma_1, omega_pow);
    HF s_sigma_2_eval = poly_eval(&s_sigma_2, omega_pow);
    HF s_sigma_3_eval = poly_eval(&s_sigma_3, omega_pow);

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
  POLY acc_x = interpolate_at_h(plonk, acc, plonk->h_len);;
  free(acc);

  // check that acc_x evaluated at omega^n is 1
  HF omega_n = hf_pow(omega, n);
  HF acc_x_eval = poly_eval(&acc_x, omega_n);
  assert(hf_equal(acc_x_eval, hf_one()));

  // create z_x polynomial
  HF b7 = rand[6], b8 = rand[7], b9 = rand[8];
  HF z_blinding_coeffs[] = {b9, b8, b7};
  POLY z_blinding_poly = poly_new(z_blinding_coeffs, 3);
  POLY z_blinded = poly_mul(&z_blinding_poly, &plonk->z_h_x);
  POLY z_x = poly_add(&z_blinded, &acc_x);

  // output of second step
  // evaluate z_x at s
  G1 z_s = srs_eval_at_s(&plonk->srs, &z_x);

  // clean up
  poly_free(&z_blinding_poly);
  poly_free(&z_blinded);

  // step 6: compute the quotient polynomial t(x)
  // compute lagrange polynomial L1(x)
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
  POLY beta_s_sigma1 = poly_scale(&s_sigma_1, beta);
  POLY a_x_beta_s_sigma1 = poly_add(&a_x, &beta_s_sigma1);
  POLY a_x_beta_s_sigma1_gamma = poly_add_hf(&a_x_beta_s_sigma1, gamma);
  POLY alpha_a_x_beta_s_sigma1_gamma = poly_scale(&a_x_beta_s_sigma1_gamma, alpha);
  POLY beta_s_sigma2 = poly_scale(&s_sigma_2, beta);
  POLY b_x_beta_s_sigma2 = poly_add(&b_x, &beta_s_sigma2);
  POLY b_x_beta_s_sigma2_gamma = poly_add_hf(&b_x_beta_s_sigma2, gamma);

  POLY beta_s_sigma3 = poly_scale(&s_sigma_3, beta);
  POLY c_x_beta_s_sigma3 = poly_add(&c_x, &beta_s_sigma3);
  POLY c_x_beta_s_sigma3_gamma = poly_add_hf(&c_x_beta_s_sigma3, gamma);

  // convert result matrix to polynomial
  HF *coeffs = (HF *)malloc(z_x.len * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < z_x.len; i++) {
    coeffs[i] = hf_mul(z_x.coeffs[i], hf_pow(omega, i));
  }
  POLY z_omega_x = poly_new(coeffs, z_x.len);
  POLY t_3_z_h = poly_mul(&alpha_a_x_beta_s_sigma1_gamma, &b_x_beta_s_sigma2_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &c_x_beta_s_sigma3_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &z_omega_x);

  poly_free(&beta_s_sigma1);
  poly_free(&a_x_beta_s_sigma1);
  poly_free(&beta_s_sigma2);
  poly_free(&b_x_beta_s_sigma2);
  poly_free(&beta_s_sigma3);
  poly_free(&c_x_beta_s_sigma3);
  free(coeffs);

  // compute t_4_z_h
  HF neg_one[1];
  neg_one[0] = hf_neg(hf_one());
  POLY _1 = poly_new(neg_one, 1);
  POLY z_x_1 = poly_add(&z_x, &_1);
  POLY alhpa_2_z_x_1 = poly_scale(&z_x_1, hf_pow(alpha, 2));
  POLY t_4_z_h = poly_mul(&alhpa_2_z_x_1, &l_1_x);

  poly_free(&_1);
  poly_free(&z_x_1);
  poly_free(&alhpa_2_z_x_1);

  POLY t_x_numer = poly_add(&t_1_z_h, &t_2_z_h);
  poly_free(&t_1_z_h);
  poly_free(&t_2_z_h);
  t_x_numer = poly_sub(&t_x_numer, &t_3_z_h);
  poly_free(&t_3_z_h);
  t_x_numer = poly_add(&t_x_numer, &t_4_z_h);
  poly_free(&t_4_z_h);

  POLY t_x;
  POLY remainder;
  poly_divide(&t_x_numer, &plonk->z_h_x, &t_x, &remainder);
  poly_free(&t_x_numer);
  if (!poly_is_zero(&remainder)) {
    fprintf(stderr, "Non-zero remainder in t(x) division\n");
    exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // step 7: split t(x) into t_lo_x, t_mid_x, t_hi_x
  size_t degree = t_x.len;
  size_t part_size = n + 2;

  POLY t_lo_x = poly_slice(&t_x, 0, part_size);
  POLY t_mid_x = poly_slice(&t_x, part_size, 2 * part_size);
  POLY t_hi_x = poly_slice(&t_x, 2 * part_size, degree);

  // step 8: evaluate t_lo_x, t_mid_x, t_hi_x at s
  G1 t_lo_s = srs_eval_at_s(&plonk->srs, &t_lo_x);
  G1 t_mid_s = srs_eval_at_s(&plonk->srs, &t_mid_x);
  G1 t_hi_s = srs_eval_at_s(&plonk->srs, &t_hi_x);

  // step 9: compute r(x) and evaluate it at z
  HF a_z = poly_eval(&a_x, z);
  HF b_z = poly_eval(&b_x, z);
  HF c_z = poly_eval(&c_x, z);
  HF s_sigma_1_z = poly_eval(&s_sigma_1, z);
  HF s_sigma_2_z = poly_eval(&s_sigma_2, z);
  HF t_z = poly_eval(&t_x, z);
  HF z_omega_z = poly_eval(&z_omega_x, z);
  poly_free(&t_x);

  // compute r_1_x
  POLY a_z_b_z_q_m_x = poly_scale(&q_m_x, hf_mul(a_z, b_z));
  POLY a_z_q_l_x = poly_scale(&q_l_x, a_z);
  POLY b_z_q_r_x = poly_scale(&q_r_x, b_z);
  POLY c_z_q_o_x = poly_scale(&q_o_x, c_z);
  POLY r_1_x = poly_add(&a_z_b_z_q_m_x, &a_z_q_l_x);
  r_1_x = poly_add(&r_1_x, &b_z_q_r_x);
  r_1_x = poly_add(&r_1_x, &c_z_q_o_x);

  // compute r_2_x
  HF a_z_beta_z_gamma = hf_add(hf_add(a_z, hf_mul(beta, z)), gamma);
  HF b_z_beta_k1_z_gamma = hf_add(hf_add(b_z, hf_mul(hf_mul(beta, k1), z)), gamma);
  HF c_z_beta_k2_z_gamma = hf_add(hf_add(c_z, hf_mul(hf_mul(beta, k2), z)), gamma);
  POLY r_2_x = poly_scale(&z_x,
                          hf_mul(
                              hf_mul(
                                  hf_mul(a_z_beta_z_gamma, b_z_beta_k1_z_gamma),
                                  c_z_beta_k2_z_gamma),
                              alpha));

  // compute r_3_x
  POLY s_sigma_3_beta_z_omega_z = poly_scale(&s_sigma_3, hf_mul(beta, z_omega_z));
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma));
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma));
  POLY r_3_x = poly_mul(&z_x, &s_sigma_3_beta_z_omega_z);
  r_3_x = poly_scale(&r_3_x,
                     hf_mul(
                         hf_mul(a_z_beta_s_sigma_1_z_gamma, b_z_beta_s_sigma_2_z_gamma),
                         alpha));

  // compute r_4_x
  POLY r_4_x = poly_scale(&z_x, hf_mul(poly_eval(&l_1_x, z), hf_pow(alpha, 2)));

  POLY r_x = poly_add(&r_1_x, &r_2_x);
  r_x = poly_add(&r_x, &r_3_x);
  r_x = poly_add(&r_x, &r_4_x);

  // compute linearization evaluation
  HF r_z = poly_eval(&r_x, z);

  // round 5
  // ---------------------------------------------------------------------------
  // we create two large polynomials that combine all the polynomials we've been
  // using so far and we output commitments to them.

  // compute opening proof polynomial w_z_x
  POLY t_mid_x_z = poly_scale(&t_mid_x, hf_pow(z, n+2));
  POLY t_hi_x_z = poly_scale(&t_hi_x, hf_pow(z, 2*n+4));
  POLY w_z_x = poly_add(&t_lo_x, &t_mid_x_z);
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
  POLY s_sigma_1_1_z_v = poly_add_hf(&s_sigma_1, hf_neg(s_sigma_1_z));
  s_sigma_1_1_z_v = poly_scale(&s_sigma_1_1_z_v, hf_pow(v, 5));
  POLY s_sigma_2_2_z_v = poly_add_hf(&s_sigma_2, hf_neg(s_sigma_2_z));
  s_sigma_2_2_z_v = poly_scale(&s_sigma_2_2_z_v, hf_pow(v, 6));
  w_z_x = poly_add(&w_z_x, &r_x_r_z_v);
  w_z_x = poly_add(&w_z_x, &a_x_a_z_v);
  w_z_x = poly_add(&w_z_x, &b_x_b_z_v);
  w_z_x = poly_add(&w_z_x, &c_x_c_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_1_1_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_2_2_z_v);

  HF coeffs_denom1[] = {hf_neg(z), hf_one()};
  POLY denom1 = poly_new(coeffs_denom1, 2);
  POLY w_z_x_quo, rem1;
  poly_divide(&w_z_x, &denom1, &w_z_x_quo, &rem1);
  assert(poly_is_zero(&rem1));

  POLY z_x_z_omega_z = poly_add_hf(&z_x, hf_neg(z_omega_z));
  HF coeffs_denom2[] = {hf_mul(hf_neg(z), omega), hf_one()};
  POLY denom2 = poly_new(coeffs_denom2, 2);
  POLY w_z_omega_x, rem2;
  poly_divide(&z_x_z_omega_z, &denom2, &w_z_omega_x, &rem2);
  assert(poly_is_zero(&rem2));

  // compute opening proof polinomials at s
  G1 w_z_s = srs_eval_at_s(&plonk->srs, &w_z_x_quo);
  G1 w_z_omega_s = srs_eval_at_s(&plonk->srs, &w_z_omega_x);

  PROOF proof;
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

  poly_free(&t_lo_x);
  poly_free(&t_mid_x);
  poly_free(&t_hi_x);
  poly_free(&z_x);
  poly_free(&l_1_x);
  poly_free(&z_omega_x);
  poly_free(&r_1_x);
  poly_free(&r_2_x);
  poly_free(&r_3_x);
  poly_free(&r_4_x);
  poly_free(&r_x);
  poly_free(&t_mid_x_z);
  poly_free(&t_hi_x_z);

  return proof;
}


#endif // PLONK_H
