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
  hf_fe alpha;
  hf_fe beta;
  hf_fe gamma;
  hf_fe z;
  hf_fe v;
} challenge;

typedef struct {
  g1_p a_s;
  g1_p b_s;
  g1_p c_s;
  g1_p z_s;
  g1_p t_lo_s;
  g1_p t_mid_s;
  g1_p t_hi_s;
  g1_p w_z_s;
  g1_p w_z_omega_s;
  hf_fe a_z;
  hf_fe b_z;
  hf_fe c_z;
  hf_fe s_sigma_1_z;
  hf_fe s_sigma_2_z;
  hf_fe r_z;
  hf_fe z_omega_z ;
} proof;

typedef struct {
  srs s;
  hf_fe *h;          // array of hf elements
  size_t h_len;
  hf_fe *k1_h;       // array of k1 * h elements
  hf_fe *k2_h;       // array of k2 * h elements
  poly z_h_x;        // polynomial z_h_x
} plonk;

plonk plonk_new(srs s, size_t n) {
  plonk p;
  p.s = s;
  p.h_len = n + 1;

  // init
  hf_fe OMEGA = hf_fe_new(OMEGA_VALUE);
  hf_fe K1 = hf_fe_new(K1_VALUE);
  hf_fe K2 = hf_fe_new(K2_VALUE);

  // compute h = [omega^0, omega^1, ..., omega^n]
  p.h = (hf_fe *)malloc(p.h_len * sizeof(hf_fe));
  if (!p.h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p.h_len; i++) {
       p.h[i] = hf_fe_pow(OMEGA, i);
  }

  // check that K1 and K2 are not in H
  for (size_t i = 0; i < p.h_len; i++) {
    if (hf_fe_equal(p.h[i], K1) || hf_fe_equal(p.h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }

  // compute h1_h: h1_h[i] = K1 * h[i]
  p.k1_h = (hf_fe *)malloc(p.h_len * sizeof(hf_fe));
  if (!p.k1_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p.h_len; i++) {
    p.k1_h[i] = hf_fe_mul(K1, p.h[i]);
  }

  // compute k2_h: k2_h[i] = K2 * h[i]
  p.k2_h = (hf_fe *)malloc(p.h_len * sizeof(hf_fe));
  if (!p.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p.h_len; i++) {
    p.k2_h[i] = hf_fe_mul(K2, p.h[i]);
  }

  // build h_pows matrix (vandermonde matrix) and compute its inverse
  matrix h_pows = matrix_zero(p.h_len, p.h_len);
  for (size_t r = 0; r < h_pows.m; r++) {
    hf_fe h_pow = hf_fe_new(1);
    for (size_t c = 0; c < h_pows.n; c++) {
      matrix_set(&h_pows, r, c, gf_from_hf(h_pow));
      h_pow = hf_fe_mul(h_pow, p.h[r]);
    }
  }

  matrix_free(&h_pows);

  // compute z_h_x = Poly::z(h)
  u8_fe h_gf = gf_from_hf(*p.h);
  p.z_h_x = poly_z(&h_gf, p.h_len);;

  return p;
}

void plonk_free(plonk *p) {
  srs_free(&p->s);
  if (p->h) {
    free(p->h);
    p->h = NULL;
  }

  if (p->k1_h) {
    free(p->k1_h);
    p->k1_h = NULL;
  }
  if (p->k2_h) {
    free(p->k2_h);
    p->k2_h = NULL;
  }
  poly_free(&p->z_h_x);
}

void copy_constraints_to_roots(const plonk *p, const copy_of *cp, size_t len, hf_fe *sigma) {
  for (size_t i = 0; i < len; i++) {
    size_t idx = cp[i].index - 1;
    switch (cp[i].type) {
    case COPYOF_A:
      sigma[i] = p->h[idx];
      break;
    case COPYOF_B:
      sigma[i] = p->k1_h[idx];
      break;
    case COPYOF_C:
      sigma[i] = p->k2_h[idx];
      break;
    default:
      fprintf(stderr, "Invalid copy_of type\n");
      exit(EXIT_FAILURE);
    }
  }
}

poly interpolate_at_h(const plonk *p, const hf_fe *values, size_t len) {
  // map roots of unity h from HF to GF
  u8_fe *h_gf = (u8_fe *)malloc(p->h_len * sizeof(u8_fe));
  if (!h_gf) {
    fprintf(stderr, "Memory allocation failed in interpolate_at_h\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p->h_len; i++) {
    h_gf[i] = gf_from_hf(p->h[i]); // p->h[i] is hf_fe
  }

  // build vandermonde matrix over GF
  matrix h_pows = matrix_zero(p->h_len, p->h_len);
  for (size_t i = 0; i < p->h_len; i++) {
    u8_fe x = h_gf[i];
    u8_fe x_pow = u8_fe_one(); // x^0
    for (size_t j = 0; j < p->h_len; j++) {
      matrix_set(&h_pows, i, j, x_pow); // set h_pows[i][j] = x_pow
      x_pow = u8_fe_mul(x_pow, x); // x_pow *= x
    }
  }

  // invert the vanderomnde matrix over GF
  matrix h_pows_inv = matrix_inv(&h_pows); // ensure matrix_inverse works over GF
  matrix_free(&h_pows);
  free(h_gf);

  // convert values to a column matrix over GF
  matrix vec = matrix_zero(len, 1);
  for (size_t i = 0; i < len; i++) {
    u8_fe tmp = gf_from_hf(values[i]); // values[i] is GF
    matrix_set(&vec, i , 0, tmp);
  }

  // multiply h_pows_inv * vec over GF
  matrix res = matrix_mul(&h_pows_inv, &vec);
  matrix_free(&h_pows_inv);
  matrix_free(&vec);

  // convert result matrix to polynomial coefficients over GF
  u8_fe *coeffs = (u8_fe *)malloc(res.m * sizeof(u8_fe));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed for coeffs\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }

  poly result = poly_new(coeffs, res.m);
  free(coeffs);
  matrix_free(&res);

  return result;
}

proof plonk_prove(
    plonk *p,
    constraints *c,
    assignments *a,
    challenge *ch,
    hf_fe rand[9]
                  ) {
  // step 1: Check that the constraints satisfy the assignments
  assert(constraints_satisfy(c, a));

  // extract challenges
  hf_fe alpha = ch->alpha;
  hf_fe beta = ch->beta;
  hf_fe gamma = ch->gamma;
  hf_fe z = ch->z;
  hf_fe v = ch->v;

  // constants
  hf_fe omega = hf_fe_new(OMEGA_VALUE);
  hf_fe k1 = hf_fe_new(K1_VALUE);
  hf_fe k2 = hf_fe_new(K2_VALUE);

  size_t n = c->num_constraints;

  // step 2: create sigmas
  hf_fe *sigma_1 = (hf_fe *)malloc(n * sizeof(hf_fe));
  hf_fe *sigma_2 = (hf_fe *)malloc(n * sizeof(hf_fe));
  hf_fe *sigma_3 = (hf_fe *)malloc(n * sizeof(hf_fe));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_prove\n");
    exit(EXIT_FAILURE);
  }

  copy_constraints_to_roots(p, c->c_a, n, sigma_1);
  copy_constraints_to_roots(p, c->c_b, n, sigma_2);
  copy_constraints_to_roots(p, c->c_c, n, sigma_3);

  // step 3: create polynomials

  // Pad a->a to length p->h_len
  hf_fe *a_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  if (!a_padded) {
    fprintf(stderr, "Memory allocation failed for a_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    a_padded[i] = f17(a->a[i].value);
  }

  poly f_a_x = interpolate_at_h(p, a_padded, p->h_len);
  free(a_padded);

  hf_fe *b_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  if (!b_padded) {
    fprintf(stderr, "Memory allocation failed for b_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    b_padded[i] = f17(a->b[i].value);
  }

  poly f_b_x = interpolate_at_h(p, b_padded, p->h_len);
  free(b_padded);

  hf_fe *c_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  if (!c_padded) {
    fprintf(stderr, "Memory allocation failed for c_padded\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    c_padded[i] = f17(a->c[i].value);
  }
  poly f_c_x = interpolate_at_h(p, c_padded, p->h_len);
  free(c_padded);

  hf_fe *q_o_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *q_m_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *q_l_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *q_r_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *q_c_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  if (!q_o_padded || !q_m_padded || !q_l_padded || !q_r_padded || !q_c_padded) {
       fprintf(stderr, "Memory allocation failed for padded constraint arrays\n");
       exit(EXIT_FAILURE);
  }
  memcpy(q_o_padded, c->q_o, n * sizeof(u8_fe));
  memcpy(q_m_padded, c->q_m, n * sizeof(u8_fe));
  memcpy(q_l_padded, c->q_l, n * sizeof(u8_fe));
  memcpy(q_r_padded, c->q_r, n * sizeof(u8_fe));
  memcpy(q_c_padded, c->q_c, n * sizeof(u8_fe));

  poly q_o_x = interpolate_at_h(p, q_o_padded, p->h_len);
  poly q_m_x = interpolate_at_h(p, q_m_padded, p->h_len);
  poly q_l_x = interpolate_at_h(p, q_l_padded, p->h_len);
  poly q_r_x = interpolate_at_h(p, q_r_padded, p->h_len);
  poly q_c_x = interpolate_at_h(p, q_c_padded, p->h_len);

  free(q_o_padded);
  free(q_m_padded);
  free(q_l_padded);
  free(q_r_padded);
  free(q_c_padded);

  // Similarly, pad sigma_1, sigma_2, sigma_3 before interpolating
  hf_fe *sigma_1_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *sigma_2_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  hf_fe *sigma_3_padded = (hf_fe *)calloc(p->h_len, sizeof(hf_fe));
  if (!sigma_1_padded || !sigma_2_padded || !sigma_3_padded) {
    fprintf(stderr, "Memory allocation failed for sigma_padded arrays\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < n; i++) {
    sigma_1_padded[i] = sigma_1[i];
    sigma_2_padded[i] = sigma_2[i];
    sigma_3_padded[i] = sigma_3[i];
  }

  poly s_sigma_1 = interpolate_at_h(p, sigma_1_padded, p->h_len);
  poly s_sigma_2 = interpolate_at_h(p, sigma_2_padded, p->h_len);
  poly s_sigma_3 = interpolate_at_h(p, sigma_3_padded, p->h_len);
  free(sigma_1_padded);
  free(sigma_2_padded);
  free(sigma_3_padded);

  // step 4: round 1 - evaluate a(x), b(x), c(x) at s
  hf_fe b1 = rand[0], b2 = rand[1], b3 = rand[2], b4 = rand[3], b5 = rand[4], b6 = rand[5];

  // a_x = (b2 + b1*x) * z_h_x + f_a_x
  hf_fe a_blinding_coeffs[] = {b2, b1};
  poly a_blinding_poly  = poly_new(a_blinding_coeffs, 2);
  poly a_x_blinded = poly_mul(&a_blinding_poly, &p->z_h_x);
  poly a_x = poly_add(&a_x_blinded, &f_a_x);

  // b_x = (b4 + b3*x) * z_h_x + f_b_x
  hf_fe b_blinding_coeffs[] = {b4, b3};
  poly b_blinding_poly  = poly_new(b_blinding_coeffs, 2);
  poly b_x_blinded = poly_mul(&b_blinding_poly, &p->z_h_x);
  poly b_x = poly_add(&b_x_blinded, &f_b_x);

  // c_x = (b6 + b5*x) * z_h_x + f_c_x
  hf_fe c_blinding_coeffs[] = {b6, b5};
  poly c_blinding_poly  = poly_new(c_blinding_coeffs, 2);
  poly c_x_blinded = poly_mul(&c_blinding_poly, &p->z_h_x);
  poly c_x = poly_add(&c_x_blinded, &f_c_x);

  // evaluate at s
  g1_p a_s = srs_eval_at_s(&p->s, &a_x);
  g1_p b_s = srs_eval_at_s(&p->s, &b_x);
  g1_p c_s = srs_eval_at_s(&p->s, &c_x);

  // cleaning up
  poly_free(&a_blinding_poly);
  poly_free(&a_x_blinded);
  poly_free(&b_blinding_poly);
  poly_free(&b_x_blinded);
  poly_free(&c_blinding_poly);
  poly_free(&c_x_blinded);

  // step 5: round 2 - evaluate the accumultor vector polynomial at s
  hf_fe b7 = rand[6], b8 = rand[7], b9 = rand[8];

  // initialize accumulator vector with 1
  u8_fe *acc = (u8_fe *)malloc(p->h_len *sizeof(u8_fe));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    exit(EXIT_FAILURE);
  }
  acc[0] = u8_fe_new(1);

  for (size_t i = 1; i <= n; i++) {
    // assignments at positoin i - 1
    hf_fe as_a = a->a[i - 1];
    hf_fe as_b = a->b[i - 1];
    hf_fe as_c = a->c[i - 1];

    // omega_pow = omega^(i-1)
    hf_fe omega_pow = hf_fe_pow(omega, i - 1);

    // denominator and numerator computations
    hf_fe denom = hf_fe_mul(
        hf_fe_mul(
            hf_fe_add(as_a, hf_fe_add(hf_fe_mul(beta, omega_pow), gamma)),
            hf_fe_add(as_b, hf_fe_add(hf_fe_mul(beta, hf_fe_mul(k1, omega_pow)), gamma))
                  ),
        u8_fe_add(as_c, u8_fe_add(u8_fe_mul(beta, u8_fe_mul(k2, omega_pow)), gamma))
			    );

    // s_sigma evaluations at omega_pow
    u8_fe s_sigma_1_eval = poly_eval(&s_sigma_1, f101(omega_pow.value));
    u8_fe s_sigma_2_eval = poly_eval(&s_sigma_2, f101(omega_pow.value));
    u8_fe s_sigma_3_eval = poly_eval(&s_sigma_3, f101(omega_pow.value));

    hf_fe numer = hf_fe_mul(
        hf_fe_mul(
            hf_fe_add(as_a, hf_fe_add(hf_fe_mul(beta, s_sigma_1_eval), gamma)),
            hf_fe_add(as_b, hf_fe_add(hf_fe_mul(beta, s_sigma_2_eval), gamma))
                  ),
        hf_fe_add(as_c, hf_fe_add(u8_fe_mul(beta, s_sigma_3_eval), gamma))
                            );

    u8_fe fraction = u8_fe_div(denom, numer);
    acc[i] = u8_fe_mul(acc[i - 1], fraction);
  }

  // pad the rest with zeros
  for (size_t i = n + 1; i < p->h_len; i++)
    acc[i] = u8_fe_new(0);

  // interpolate accumualtor vector
  poly acc_x = interpolate_at_h(p, acc, p->h_len);
  free(acc);

  // check that acc_x evaluated at omega^n is 1
  u8_fe omega_n = u8_fe_pow(omega, n);
  printf("omega_n = %u\n", omega_n.value);
  u8_fe acc_x_eval = poly_eval(&acc_x, omega_n);
  printf("acc_x_eval = %u\n", acc_x_eval.value);
  printf("Evaluating acc_x at x = %u\n", omega_n.value);
  printf("acc_x coefficients:\n");
  for (size_t i = 0; i < acc_x.len; i++) {
    printf("Coefficient of x^%zu: %u\n", i, acc_x.coeffs[i].value);
  }
  printf("acc values:\n");
  for (size_t i = 0; i < p->h_len; i++) {
    printf("acc[%zu] = %u\n", i, acc[i].value);
  }

  assert(u8_fe_equal(acc_x_eval, u8_fe_new(1)));

  // create z_x polynomial
  u8_fe z_blinding_coeffs[] = {b9, b8, b7};
  poly z_blinding_poly = poly_new(z_blinding_coeffs, 3);
  poly z_blinded = poly_mul(&z_blinding_poly, &p->z_h_x);
  poly z_x = poly_add(&z_blinded, &acc_x);

  // evaluate z_x at s
  g1_p z_s = srs_eval_at_s(&p->s, &z_x);

  // clean up
  poly_free(&z_blinding_poly);
  poly_free(&z_blinded);

  // step 6: Compute the quotient polynomial t(x)
  // compute Lagrange polynomial L1(x)
  u8_fe *lagrange_vector = (u8_fe *)calloc(p->h_len, sizeof(u8_fe));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed for lagrange_vector\n");
    exit(EXIT_FAILURE);
  }
  lagrange_vector[0] = u8_fe_new(1); // L1(omega^0) = 1, rest are 0
  poly l_1_x = interpolate_at_h(p, lagrange_vector, p->h_len);
  free(lagrange_vector);

  // compute p_i_x (public input polynomial)
  // assuming no public inputs for simplicity
  poly p_i_x = poly_zero();

  // compute helper polynomials
  poly a_x_b_x = poly_mul(&a_x, &b_x);
  poly a_x_b_x_q_m_x = poly_mul(&a_x_b_x, &q_m_x);
  poly_free(&a_x_b_x);

  poly a_x_q_l_x = poly_mul(&a_x, &q_l_x);
  poly b_x_q_r_x = poly_mul(&b_x, &q_r_x);
  poly c_x_q_o_x = poly_mul(&c_x, &q_o_x);

  poly sum1 = poly_add(&a_x_b_x_q_m_x, &a_x_q_l_x);
  poly_free(&a_x_b_x_q_m_x);
  poly_free(&a_x_q_l_x);

  poly sum2 = poly_add(&b_x_q_r_x, &c_x_q_o_x);
  poly_free(&b_x_q_r_x);
  poly_free(&c_x_q_o_x);

  poly t_1_z_h = poly_add(&sum1, &sum2);
  poly_free(&sum1);
  poly_free(&sum2);

  t_1_z_h = poly_add(&t_1_z_h, &p_i_x);
  poly_free(&p_i_x);
  t_1_z_h = poly_add(&t_1_z_h, &q_c_x);

  // compute t_2_z_h
  poly tmp = poly_new(&gamma, 1);
  poly a_x_beta = poly_add(&a_x, &tmp);
  tmp = poly_new(&omega, 1);
  tmp = poly_scale(&tmp, beta);
  a_x_beta = poly_add(&a_x_beta, &tmp);
  tmp = poly_new(&gamma, 1);
  poly b_x_beta_k1 = poly_add(&b_x, &tmp);
  u8_fe beta_k1 = u8_fe_mul(beta, k1);
  poly beta_k1_1 = poly_new(&beta_k1, 1);
  poly beta_k1_1_1 = poly_scale(&beta_k1_1, u8_fe_new(1));
  b_x_beta_k1 = poly_add(&b_x_beta_k1, &beta_k1_1_1);
  poly c_x_beta_k2 = poly_add(&c_x, &tmp);
  u8_fe beta_k2 = u8_fe_mul(beta, k2);
  poly beta_k2_1 = poly_new(&beta_k2, 1);
  poly beta_k2_1_1 = poly_scale(&beta_k2_1, u8_fe_new(1));
  c_x_beta_k2 = poly_add(&c_x_beta_k2, &beta_k2_1_1);
  poly temp = poly_mul(&a_x_beta, &b_x_beta_k1);
  poly t_2_z_h = poly_mul(&temp, &c_x_beta_k2);
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
  poly s_sigma_1_beta = poly_scale(&s_sigma_1, beta);
  tmp = poly_new(&gamma, 1);
  s_sigma_1_beta = poly_add(&s_sigma_1_beta, &tmp);
  poly s_sigma_2_beta = poly_scale(&s_sigma_2, beta);
  s_sigma_2_beta = poly_add(&s_sigma_2_beta, &tmp);
  poly s_sigma_3_beta = poly_scale(&s_sigma_3, beta);
  s_sigma_3_beta = poly_add(&s_sigma_3_beta, &tmp);
  poly_free(&tmp);

  poly temp1 = poly_add(&a_x, &s_sigma_1_beta);
  poly temp2 = poly_add(&b_x, &s_sigma_2_beta);
  poly temp3 = poly_add(&c_x, &s_sigma_3_beta);
  poly_free(&s_sigma_1_beta);
  poly_free(&s_sigma_2_beta);
  poly_free(&s_sigma_3_beta);

  poly temp4 = poly_mul(&temp1, &temp2);
  poly temp5 = poly_mul(&temp4, &temp3);
  poly_free(&temp1);
  poly_free(&temp2);
  poly_free(&temp3);
  poly_free(&temp4);

  poly z_omega_x = poly_shift(&z_x, 1); // Shift z_x by multiplying x^1
  poly t_3_z_h = poly_mul(&temp5, &z_omega_x);
  poly_free(&temp5);
  poly_free(&z_omega_x);
  t_3_z_h = poly_scale(&t_3_z_h, alpha);

  // compute t_4_z_h
  u8_fe one = u8_fe_one();
  poly one_1 = poly_new(&one, 1);
  poly one_1_neg = poly_negate(&one_1);
  poly z_x_minus_one = poly_add(&z_x, &one_1_neg);
  poly t_4_z_h = poly_scale(&z_x_minus_one, u8_fe_pow(alpha, 2));
  poly_free(&one_1);
  poly_free(&one_1_neg);
  poly_free(&z_x_minus_one);
  t_4_z_h = poly_mul(&t_4_z_h, &l_1_x);

  // compute t(x)
  poly t_x_numerator = poly_add(&t_1_z_h, &t_2_z_h);
  poly_free(&t_1_z_h);
  poly_free(&t_2_z_h);
  t_x_numerator = poly_sub(&t_x_numerator, &t_3_z_h);
  poly_free(&t_3_z_h);
  t_x_numerator = poly_add(&t_x_numerator, &t_4_z_h);
  poly_free(&t_4_z_h);

  // divide t_x_numerator by z_h_x to get t_x
  poly t_x;
  poly remainder;
  poly_divide(&t_x_numerator, &p->z_h_x, &t_x, &remainder);
  poly_free(&t_x_numerator);
  if (!poly_is_zero(&remainder)) {
    fprintf(stderr, "Non-zero remainder in t(x) division\n");
    exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // evaluate t_x at z to get t_z
  u8_fe t_z = poly_eval(&t_x, z);

  // step 7: Split t(x) into t_lo_x, t_mid_x, t_hi_x
  size_t degree = t_x.len;
  size_t n_plus_2 = n + 2;
  size_t part_size = n_plus_2;
  size_t total_parts = 3;

  // ensure t_x has enough coefficients
  if (degree < 3 * part_size) {
    // pad t_x with zeros if necessary
    u8_fe *padded_coeffs = (u8_fe *)calloc(3 * part_size, sizeof(u8_fe));
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

  poly t_lo_x = poly_slice(&t_x, 0, part_size);
  poly t_mid_x = poly_slice(&t_x, part_size, 2 * part_size);
  poly t_hi_x = poly_slice(&t_x, 2 * part_size, degree);
  poly_free(&t_x);

  // step 8: Evaluate t_lo_x, t_mid_x, t_hi_x at s
  g1_p t_lo_s = srs_eval_at_s(&p->s, &t_lo_x);
  g1_p t_mid_s = srs_eval_at_s(&p->s, &t_mid_x);
  g1_p t_hi_s = srs_eval_at_s(&p->s, &t_hi_x);

  // clean up t polynomials
  poly_free(&t_lo_x);
  poly_free(&t_mid_x);
  poly_free(&t_hi_x);

  // Step 9: Compute r(x) and evaluate it at z
  // Evaluate polynomials at z
  u8_fe a_z = poly_eval(&a_x, z);
  u8_fe b_z = poly_eval(&b_x, z);
  u8_fe c_z = poly_eval(&c_x, z);
  u8_fe s_sigma_1_z = poly_eval(&s_sigma_1, z);
  u8_fe s_sigma_2_z = poly_eval(&s_sigma_2, z);
  u8_fe z_omega_z = poly_eval(&z_x, u8_fe_mul(z, omega));

  // Compute r_1_x
  poly q_m_x_a_z_b_z = poly_scale(&q_m_x, u8_fe_mul(a_z, b_z));
  poly q_l_x_a_z = poly_scale(&q_l_x, a_z);
  poly q_r_x_b_z = poly_scale(&q_r_x, b_z);
  poly q_o_x_c_z = poly_scale(&q_o_x, c_z);
  poly sum_r1 = poly_add(&q_m_x_a_z_b_z, &q_l_x_a_z);
  poly_free(&q_m_x_a_z_b_z);
  poly_free(&q_l_x_a_z);
  sum_r1 = poly_add(&sum_r1, &q_r_x_b_z);
  poly_free(&q_r_x_b_z);
  sum_r1 = poly_add(&sum_r1, &q_o_x_c_z);
  poly_free(&q_o_x_c_z);
  poly r_1_x = poly_add(&sum_r1, &q_c_x);
  poly_free(&sum_r1);

  // Compute r_2_x
  u8_fe beta_z = u8_fe_mul(beta, z);
  u8_fe gamma_beta_z = u8_fe_add(gamma, beta_z);
  u8_fe beta_k1_z = u8_fe_mul(beta, u8_fe_mul(k1, z));
  u8_fe gamma_beta_k1_z = u8_fe_add(gamma, beta_k1_z);
  u8_fe beta_k2_z = u8_fe_mul(beta, u8_fe_mul(k2, z));
  u8_fe gamma_beta_k2_z = u8_fe_add(gamma, beta_k2_z);

  u8_fe term1 = u8_fe_add(a_z, gamma_beta_z);
  u8_fe term2 = u8_fe_add(b_z, gamma_beta_k1_z);
  u8_fe term3 = u8_fe_add(c_z, gamma_beta_k2_z);

  u8_fe numerator = u8_fe_mul(u8_fe_mul(term1, term2), term3);
  poly r_2_x = poly_scale(&z_x, u8_fe_mul(numerator, alpha));

  // Compute r_3_x
  u8_fe beta_s_sigma_1_z = u8_fe_mul(beta, s_sigma_1_z);
  u8_fe gamma_beta_s_sigma_1_z = u8_fe_add(gamma, beta_s_sigma_1_z);
  u8_fe beta_s_sigma_2_z = u8_fe_mul(beta, s_sigma_2_z);
  u8_fe gamma_beta_s_sigma_2_z = u8_fe_add(gamma, beta_s_sigma_2_z);

  u8_fe term4 = u8_fe_add(a_z, gamma_beta_s_sigma_1_z);
  u8_fe term5 = u8_fe_add(b_z, gamma_beta_s_sigma_2_z);
  u8_fe temp_alpha = u8_fe_mul(alpha, alpha);

  u8_fe s_sigma_3_z = poly_eval(&s_sigma_3, z);
  u8_fe beta_s_sigma_3_z = u8_fe_mul(beta, s_sigma_3_z);
  u8_fe gamma_beta_s_sigma_3_z = u8_fe_add(gamma, beta_s_sigma_3_z);

  u8_fe term6 = u8_fe_mul(z_omega_z, beta_s_sigma_3_z);

  u8_fe numerator_r3 = u8_fe_mul(u8_fe_mul(term4, term5), term6);
  poly r_3_x = poly_scale(&z_x, u8_fe_mul(numerator_r3, alpha));

  // Compute r_4_x
  u8_fe l_1_z = poly_eval(&l_1_x, z);
  poly r_4_x = poly_scale(&z_x, u8_fe_mul(u8_fe_mul(temp_alpha, l_1_z), u8_fe_new(1)));

  // Compute r_x
  poly sum_r = poly_add(&r_1_x, &r_2_x);
  poly_free(&r_1_x);
  poly_free(&r_2_x);
  sum_r = poly_add(&sum_r, &r_3_x);
  poly_free(&r_3_x);
  sum_r = poly_add(&sum_r, &r_4_x);
  poly_free(&r_4_x);

  poly r_x = sum_r;

  // Evaluate r_x at z to get r_z
  u8_fe r_z = poly_eval(&r_x, z);

  // Step 10: Compute opening prf w_z_x and w_z_omega_x
  // Compute t_lo_x + t_mid_x * z^{n+2} + t_hi_x * z^{2n+4}
  u8_fe z_pow_n_plus_2 = u8_fe_pow(z, n + 2);
  u8_fe z_pow_2n_plus_4 = u8_fe_pow(z, 2 * n + 4);

  poly t_mid_x_z = poly_scale(&t_mid_x, z_pow_n_plus_2);
  poly t_hi_x_z = poly_scale(&t_hi_x, z_pow_2n_plus_4);

  poly t_sum = poly_add(&t_lo_x, &t_mid_x_z);
  poly_free(&t_mid_x_z);
  t_sum = poly_add(&t_sum, &t_hi_x_z);
  poly_free(&t_hi_x_z);

  // Compute w_z_x numerator
  poly t_z_poly = poly_new(&t_z, 1); // Polynomial representing the constant t_z
  poly w_z_x_num = poly_sub(&t_sum, &t_z_poly);
  poly_free(&t_sum);
  poly_free(&t_z_poly);

  // Compute (r_x - r_z) * v
  poly r_z_1 = poly_new(&r_z, 1);
  poly r_x_minus_r_z = poly_sub(&r_x, &r_z_1);
  poly r_x_minus_r_z_v = poly_scale(&r_x_minus_r_z, v);
  poly_free(&r_z_1);
  poly_free(&r_x_minus_r_z);

  // Compute (a_x - a_z) * v^2
  poly a_z_1 = poly_new(&a_z, 1);
  poly a_x_minus_a_z = poly_sub(&a_x, &a_z_1);
  u8_fe v_squared = u8_fe_mul(v, v);
  poly a_x_minus_a_z_v2 = poly_scale(&a_x_minus_a_z, v_squared);
  poly_free(&a_z_1);
  poly_free(&a_x_minus_a_z);

  // Similarly compute (b_x - b_z) * v^3 and (c_x - c_z) * v^4
  u8_fe v_cubed = u8_fe_mul(v_squared, v);
  poly b_z_1 = poly_new(&b_z, 1);
  poly b_x_minus_b_z = poly_sub(&b_x, &b_z_1);
  poly b_x_minus_b_z_v3 = poly_scale(&b_x_minus_b_z, v_cubed);
  poly_free(&b_z_1);
  poly_free(&b_x_minus_b_z);

  u8_fe v_fourth = u8_fe_mul(v_cubed, v);
  poly c_z_1 = poly_new(&c_z, 1);
  poly c_x_minus_c_z = poly_sub(&c_x, &c_z_1);
  poly c_x_minus_c_z_v4 = poly_scale(&c_x_minus_c_z, v_fourth);
  poly_free(&c_z_1);
  poly_free(&c_x_minus_c_z);

  // Similarly compute terms for s_sigma_1 and s_sigma_2
  u8_fe v_fifth = u8_fe_mul(v_fourth, v);
  poly s_sigma_1_z_1 = poly_new(&s_sigma_1_z, 1);
  poly s_sigma_1_x_minus_s_sigma_1_z = poly_sub(&s_sigma_1, &s_sigma_1_z_1);
  poly s_sigma_1_x_minus_s_sigma_1_z_v5 = poly_scale(&s_sigma_1_x_minus_s_sigma_1_z, v_fifth);
  poly_free(&s_sigma_1_z_1);
  poly_free(&s_sigma_1_x_minus_s_sigma_1_z);

  u8_fe v_sixth = u8_fe_mul(v_fifth, v);
  poly s_sigma_2_z_1 = poly_new(&s_sigma_2_z, 1);
  poly s_sigma_2_x_minus_s_sigma_2_z = poly_sub(&s_sigma_2, &s_sigma_2_z_1);
  poly s_sigma_2_x_minus_s_sigma_2_z_v6 = poly_scale(&s_sigma_2_x_minus_s_sigma_2_z, v_sixth);
  poly_free(&s_sigma_2_z_1);
  poly_free(&s_sigma_2_x_minus_s_sigma_2_z);

  // Sum all terms to compute w_z_x numerator
  poly w_z_x_num_total = poly_add(&w_z_x_num, &r_x_minus_r_z_v);
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
  poly x_minus_z = poly_new((u8_fe[]){u8_fe_neg(z), u8_fe_new(1)}, 2);
  poly w_z_x;
  poly rem;
  poly_divide(&w_z_x_num_total, &x_minus_z, &w_z_x, &rem);
  poly_free(&w_z_x_num_total);
  poly_free(&x_minus_z);
  if (!poly_is_zero(&remainder)) {
       fprintf(stderr, "Non-zero remainder in w_z_x division\n");
       exit(EXIT_FAILURE);
  }
  poly_free(&remainder);

  // Compute w_z_s = srs_eval_at_s(w_z_x)
  g1_p w_z_s = srs_eval_at_s(&p->s, &w_z_x);

  // Compute w_z_omega_x
  poly z_omega_z_1 = poly_new(&z_omega_z, 1);
  poly z_x_minus_z_omega_z = poly_sub(&z_x, &z_omega_z_1);
  u8_fe z_omega = u8_fe_mul(z, omega);
  poly x_minus_z_omega = poly_new((u8_fe[]){u8_fe_neg(z_omega), u8_fe_new(1)}, 2);
  poly w_z_omega_x;
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
  g1_p w_z_omega_s = srs_eval_at_s(&p->s, &w_z_omega_x);

  // Clean up polynomials
  poly_free(&w_z_x);
  poly_free(&w_z_omega_x);

  // Step 11: Assemble the Proof
  proof prf = {
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
    plonk *p,
    constraints *cons,
    proof *prf,
    challenge *ch,
    u8_fe rand[1]
                  ) {
  // unpack prf
  g1_p a_s = prf->a_s;
  g1_p b_s = prf->b_s;
  g1_p c_s = prf->c_s;
  g1_p z_s = prf->z_s;
  g1_p t_lo_s = prf->t_lo_s;
  g1_p t_mid_s = prf->t_mid_s;
  g1_p t_hi_s = prf->t_hi_s;
  g1_p w_z_s = prf->w_z_s;
  g1_p w_z_omega_s = prf->w_z_omega_s;
  u8_fe a_z = prf->a_z;
  u8_fe b_z = prf->b_z;
  u8_fe c_z = prf->c_z;
  u8_fe s_sigma_1_z = prf->s_sigma_1_z;
  u8_fe s_sigma_2_z = prf->s_sigma_2_z;
  u8_fe r_z = prf->r_z;
  u8_fe z_omega_z = prf->z_omega_z;

  // unpack challenge
  u8_fe alpha = ch->alpha;
  u8_fe beta = ch->beta;
  u8_fe gamma = ch->gamma;
  u8_fe z = ch->z;
  u8_fe v = ch->v;

  // constants
  u8_fe omega = u8_fe_new(OMEGA_VALUE);
  u8_fe k1 = u8_fe_new(K1_VALUE);
  u8_fe k2 = u8_fe_new(K2_VALUE);

  // step 1: Validate proof points in G1
  if (!g1_p_is_on_curve(&a_s) ||
      !g1_p_is_on_curve(&b_s) ||
      !g1_p_is_on_curve(&c_s) ||
      !g1_p_is_on_curve(&z_s) ||
      !g1_p_is_on_curve(&t_lo_s) ||
      !g1_p_is_on_curve(&t_mid_s) ||
      !g1_p_is_on_curve(&t_hi_s) ||
      !g1_p_is_on_curve(&w_z_s) ||
      !g1_p_is_on_curve(&w_z_omega_s)) {
    return false;
  }

  // step 2: Validate proof fields in HF
  if (!u8_fe_in_field(a_z) ||
      !u8_fe_in_field(b_z) ||
      !u8_fe_in_field(c_z) ||
      !u8_fe_in_field(s_sigma_1_z) ||
      !u8_fe_in_field(s_sigma_2_z) ||
      !u8_fe_in_field(r_z) ||
      !u8_fe_in_field(z_omega_z)) {
    return false;
  }

  // step 3: No public inputs, nothing to do

  // step 4: Evaluate z_h at z
  u8_fe z_h_z = poly_eval(&p->z_h_x, z);

  // step 5: Evaluate lagrange polynomial L1 at z
  u8_fe *lagrange_vector = (u8_fe *)calloc(p->h_len, sizeof(u8_fe));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    exit(EXIT_FAILURE);
  }

  lagrange_vector[0] = u8_fe_new(1); // L1(omega^0) = 1
  poly l_1_x = interpolate_at_h(p, lagrange_vector, p->h_len);
  free(lagrange_vector);
  u8_fe l_1_z = poly_eval(&l_1_x, z);
  poly_free(&l_1_x);

  // step 6: No public inputs, nothing to do
  u8_fe p_i_z = u8_fe_new(0);

  // step 7: Compute quotient polynomial evaluation
  u8_fe a_z_beta_s_sigma_1_z_gamma = u8_fe_add(u8_fe_add(u8_fe_mul(beta, s_sigma_1_z), gamma), a_z);
  u8_fe b_z_beta_s_sigma_2_z_gamma = u8_fe_add(u8_fe_add(u8_fe_mul(beta, s_sigma_2_z), gamma), b_z);
  u8_fe c_z_gamma = u8_fe_add(c_z, gamma);
  u8_fe l_1_z_alpha_2 = u8_fe_mul(l_1_z, u8_fe_pow(alpha, 2));

  u8_fe numerator = u8_fe_sub(u8_fe_add(r_z, p_i_z),
                              u8_fe_add(
                                  u8_fe_mul(
                                      u8_fe_mul(
                                          u8_fe_mul(
                                              a_z_beta_s_sigma_1_z_gamma,
                                              b_z_beta_s_sigma_2_z_gamma),
                                          c_z_gamma),
                                      z_omega_z),
                                  l_1_z_alpha_2));

  u8_fe t_z = u8_fe_div(numerator, z_h_z);

  // step 8: Compute the first part of batched polynomial commitment
  // compute sigmas
  size_t n = cons->num_constraints;

  u8_fe *sigma_1 = (u8_fe *)malloc(n * sizeof(u8_fe));
  u8_fe *sigma_2 = (u8_fe *)malloc(n * sizeof(u8_fe));
  u8_fe *sigma_3 = (u8_fe *)malloc(n * sizeof(u8_fe));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    exit(EXIT_FAILURE);
  }

  copy_constraints_to_roots(p, cons->c_a, n, sigma_1);
  copy_constraints_to_roots(p, cons->c_b, n, sigma_2);
  copy_constraints_to_roots(p, cons->c_c, n, sigma_3);

  // interpolate and evaluate polynomials at s
  poly q_m_x = interpolate_at_h(p, cons->q_m, n);
  g1_p q_m_s = srs_eval_at_s(&p->s, &q_m_x);
  poly_free(&q_m_x);

  poly q_l_x = interpolate_at_h(p, cons->q_l, n);
  g1_p q_l_s = srs_eval_at_s(&p->s, &q_l_x);
  poly_free(&q_l_x);

  poly q_r_x = interpolate_at_h(p, cons->q_r, n);
  g1_p q_r_s = srs_eval_at_s(&p->s, &q_r_x);
  poly_free(&q_r_x);

  poly q_o_x = interpolate_at_h(p, cons->q_o, n);
  g1_p q_o_s = srs_eval_at_s(&p->s, &q_o_x);
  poly_free(&q_o_x);

  poly q_c_x = interpolate_at_h(p, cons->q_c, n);
  g1_p q_c_s = srs_eval_at_s(&p->s, &q_c_x);
  poly_free(&q_c_x);

  poly s_sigma_1_x = interpolate_at_h(p, sigma_1, n);
  g1_p sigma_1_s = srs_eval_at_s(&p->s, &s_sigma_1_x);
  poly_free(&s_sigma_1_x);

  poly s_sigma_2_x = interpolate_at_h(p, sigma_2, n);
  g1_p sigma_2_s = srs_eval_at_s(&p->s, &s_sigma_2_x);
  poly_free(&s_sigma_2_x);

  poly s_sigma_3_x = interpolate_at_h(p, sigma_3, n);
  g1_p sigma_3_s = srs_eval_at_s(&p->s, &s_sigma_3_x);
  poly_free(&s_sigma_3_x);

  free(sigma_1);
  free(sigma_2);
  free(sigma_3);

  // Random scalar u
  u8_fe u = rand[0];

  // Compute d_1_s
  g1_p tmp_lhs = g1_p_mul(&q_r_s, u8_fe_mul(b_z, v).value);
  g1_p tmp_rhs = g1_p_mul(&q_o_s, u8_fe_mul(c_z, v).value);
  g1_p tmp_lhs2 = g1_p_mul(&q_l_s, u8_fe_mul(a_z, v).value);
  g1_p tmp_rhs2 = g1_p_add(&tmp_lhs, &tmp_rhs);
  g1_p tmp_lhs3 = g1_p_mul(&q_m_s, u8_fe_mul(u8_fe_mul(a_z, b_z), v).value);
  g1_p tmp_rhs3 = g1_p_add(&tmp_lhs2, &tmp_rhs2);
  g1_p d_1_s = g1_p_add(&tmp_lhs3,&tmp_rhs3);
  g1_p tmp_rhs4 = g1_p_mul(&q_c_s, v.value);
  d_1_s = g1_p_add(&d_1_s, &tmp_rhs4);

  // Compute d_2_s
  u8_fe term_d2 = u8_fe_add(
      u8_fe_mul(
          u8_fe_mul(
              u8_fe_mul(
                  u8_fe_add(a_z, u8_fe_add(u8_fe_mul(beta, z), gamma)),
                  u8_fe_add(b_z, u8_fe_add(u8_fe_mul(beta, u8_fe_mul(k1, z)), gamma))
                        ),
              u8_fe_add(c_z, u8_fe_add(u8_fe_mul(beta, u8_fe_mul(k2, z)), gamma))
                    ),
          u8_fe_mul(alpha, v)
                ),
      u8_fe_mul(l_1_z, u8_fe_mul(u8_fe_pow(alpha, 2), v))
                            );
  term_d2 = u8_fe_add(term_d2, u);
  g1_p d_2_s = g1_p_mul(&z_s, term_d2.value);

  // Compute d_3_s
  u8_fe term_d3 = u8_fe_mul(
      u8_fe_mul(
          u8_fe_mul(
              u8_fe_add(a_z, u8_fe_add(u8_fe_mul(beta, s_sigma_1_z), gamma)),
              u8_fe_add(b_z, u8_fe_add(u8_fe_mul(beta, s_sigma_2_z), gamma))
                    ),
          u8_fe_mul(alpha, v)
                ),
      u8_fe_mul(beta, z_omega_z)
                            );
  g1_p d_3_s = g1_p_mul(&sigma_3_s, term_d3.value);

  // Compute d_s = d_1_s + d_2_s - d_3_s
  g1_p d_s = g1_p_add(&d_1_s, &d_2_s);
  g1_p d_3_s_neg = g1_p_neg(&d_3_s);
  d_s = g1_p_add(&d_s, &d_3_s_neg);

  // Step 9: Compute f_s
  g1_p t_lo_s_scaled = t_lo_s;
  g1_p t_mid_s_scaled = g1_p_mul(&t_mid_s, u8_fe_pow(z, n + 2).value);
  g1_p t_hi_s_scaled = g1_p_mul(&t_hi_s, u8_fe_pow(z, 2 * n + 4).value);
  g1_p t_lo_mid_s_scaled = g1_p_add(&t_lo_s_scaled, &t_mid_s_scaled);
  g1_p t_hi_s_scaled_d_s = g1_p_add(&t_hi_s_scaled, &d_s);
  g1_p a_s_pow_v_2 = g1_p_mul(&a_s, u8_fe_pow(v, 2).value);
  g1_p b_s_pow_v_3 = g1_p_mul(&b_s, u8_fe_pow(v, 3).value);
  g1_p c_s_pow_v_4 = g1_p_mul(&c_s, u8_fe_pow(v, 4).value);
  g1_p sigma_1_s_pow_v_5 = g1_p_mul(&sigma_1_s, u8_fe_pow(v, 5).value);
  g1_p sigma_2_s_pow_v_6 = g1_p_mul(&sigma_2_s, u8_fe_pow(v, 6).value);
  g1_p temp_sum = g1_p_add(&sigma_1_s_pow_v_5, &sigma_2_s_pow_v_6);
  temp_sum = g1_p_add(&c_s_pow_v_4, &temp_sum);
  temp_sum = g1_p_add(&b_s_pow_v_3, &temp_sum);
  temp_sum = g1_p_add(&a_s_pow_v_2, &temp_sum);

  g1_p f_s = g1_p_add(&t_lo_mid_s_scaled, &t_hi_s_scaled_d_s);
  f_s = g1_p_add(&f_s, &temp_sum);

  // Step 10: Compute e_s
  poly one_poly = poly_new(&(u8_fe){ .value = 1 }, 1);
  g1_p g1_s = srs_eval_at_s(&p->s, &one_poly);
  poly_free(&one_poly);

  u8_fe e_scalar = u8_fe_add(
      t_z,
      u8_fe_add(
          u8_fe_mul(v, r_z),
          u8_fe_add(
              u8_fe_mul(u8_fe_pow(v, 2), a_z),
              u8_fe_add(
                  u8_fe_mul(u8_fe_pow(v, 3), b_z),
                  u8_fe_add(
                      u8_fe_mul(u8_fe_pow(v, 4), c_z),
                      u8_fe_add(
                          u8_fe_mul(u8_fe_pow(v, 5), s_sigma_1_z),
                          u8_fe_add(
                              u8_fe_mul(u8_fe_pow(v, 6), s_sigma_2_z),
                              u8_fe_mul(u, z_omega_z)
                                    )
                                )
                            )
                        )
                    )
                )
                             );
  g1_p e_s = g1_p_mul(&g1_s, e_scalar.value);

  // step 11: Batch validate all equations
  g1_p w_z_omega_s_u = g1_p_mul(&w_z_omega_s, u.value);
  g1_p e_1_q1 = g1_p_add(&w_z_s, &w_z_omega_s_u);
  g2_p e_1_q2 = p->s.g2_s;

  g1_p temp1 = g1_p_mul(&w_z_s, z.value);
  g1_p temp2 = g1_p_mul(&w_z_omega_s, u8_fe_mul(u, u8_fe_mul(z, omega)).value);
  g1_p temp3 = g1_p_add(&temp1, &temp2);

  // TODO: check it's validity for g1_p_sub
  g1_p e_s_neg = g1_p_neg(&e_s);
  g1_p temp4 = g1_p_add(&f_s, &e_s_neg);

  g1_p e_2_q1 = g1_p_add(&temp3, &temp4);
  g2_p e_2_q2 = p->s.g2_1;

  // compute pairings
  gtp e_1 = pairing(&e_1_q1, &e_1_q2);
  gtp e_2 = pairing(&e_2_q1, &e_2_q2);

  // check if e_1 == e_2
  if (gtp_equal(&e_1, &e_2)) {
    return true;
  } else {
    return false;
  }
}

#endif // PLONK_H
