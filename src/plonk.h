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

#define OMEGA_VALUE 3 // example value; should be a primitive roor of unity in the field
#define K1_VALUE    2 // should not be in H
#define K2_VALUE    4 // should not be in H or K1 * H

typedef struct {
  u8_fe alpha;
  u8_fe beta;
  u8_fe gamma;
  u8_fe z;
  u8_fe v;
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
  u8_fe a_z;
  u8_fe b_z;
  u8_fe c_z;
  u8_fe s_sigma_1_z;
  u8_fe s_sigma_2_z;
  u8_fe r_z;
  u8_fe z_omega_z ;
} proof;

typedef struct {
  srs s;
  u8_fe *h;          // array of hf elements
  size_t h_len;
  matrix h_pows_inv; // inverse of h powers matrix
  u8_fe *k1_h;       // array of k1 * h elements
  u8_fe *k2_h;       // array of k2 * h elements
  poly z_h_x;        // polynomial z_h_x
} plonk;

plonk plonk_new(srs s, size_t omega_pows) {
  plonk p;
  p.s = s;
  p.h_len = omega_pows;

  // init
  u8_fe OMEGA = u8_fe_new(OMEGA_VALUE);
  u8_fe K1 = u8_fe_new(K1_VALUE);
  u8_fe K2 = u8_fe_new(K2_VALUE);

  p.h = (u8_fe *)malloc(p.h_len * sizeof(u8_fe));
  if (!p.h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  u8_fe omega_pow = u8_fe_new(1);
  for (size_t i = 0; i < p.h_len; i++) {
    p.h[i] = omega_pow;
    omega_pow = u8_fe_mul(omega_pow, OMEGA);
  }

  // check that K1 and K2 are not in H
  for (size_t i = 0; i < p.h_len; i++) {
    if (u8_fe_equal(p.h[i], K1) || u8_fe_equal(p.h[i], K2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }

  // compute h1_h: h1_h[i] = K1 * h[i]
  p.k1_h = (u8_fe *)malloc(p.h_len * sizeof(u8_fe));
  if (!p.k1_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p.h_len; i++) {
    p.k1_h[i] = u8_fe_mul(K1, p.h[i]);
  }

  // compute k2_h: k2_h[i] = K2 * h[i]
  p.k2_h = (u8_fe *)malloc(p.h_len * sizeof(u8_fe));
  if (!p.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < p.h_len; i++) {
    p.k2_h[i] = u8_fe_mul(K2, p.h[i]);
  }

  // build h_pows matrix (vandermonde matrix) and compute its inverse
  matrix h_pows = matrix_zero(p.h_len, p.h_len);
  for (size_t r = 0; r < h_pows.m; r++) {
    u8_fe h_pow = u8_fe_new(1);
    for (size_t c = 0; c < h_pows.n; c++) {
      matrix_set(&h_pows, r, c, h_pow);
      h_pow = u8_fe_mul(h_pow, p.h[r]);
    }
  }

  // compute invers eof h_pows
  p.h_pows_inv = matrix_inv(&h_pows);
  matrix_free(&h_pows);

  // compute z_h_x = Poly::z(h)
  p.z_h_x = poly_z(p.h, p.h_len);;

  return p;
}

void plonk_free(plonk *p) {
  srs_free(&p->s);
  if (p->h) {
    free(p->h);
    p->h = NULL;
  }
  matrix_free(&p->h_pows_inv);
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

void copy_constraints_to_roots(const plonk *p, const copy_of *cp, size_t len, u8_fe *sigma) {
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

poly interpolate_at_h(const plonk *p, const u8_fe *values, size_t len) {
  // convert values to a column matrix
  matrix vec = matrix_zero(len, 1);
  for (size_t i = 0; i < len; i++) {
    matrix_set(&vec, i , 0, values[i]);
  }

  // multiply h_pows_inv * vec
  matrix res = matrix_mul(&p->h_pows_inv, &vec);

  // convert result matrix to polynomial
  u8_fe *coeffs = (u8_fe *)malloc(res.m * sizeof(u8_fe));
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }

  poly result = poly_new(coeffs, res.m);

  matrix_free(&vec);
  matrix_free(&res);
  free(coeffs);

  return result;
}

proof plonk_prove(
    plonk *p,
    constraints *c,
    assignments *a,
    challenge *ch,
    u8_fe rand[9]
                  ) {
  // step 1: Check that the constraints satisfy the assignments
  assert(constraints_satisfy(c, a));

  // extract challenges
  u8_fe alpha = ch->alpha;
  u8_fe beta = ch->beta;
  u8_fe gamma = ch->gamma;
  u8_fe z = ch->z;
  u8_fe v = ch->v;

  // constants
  u8_fe omega = u8_fe_new(OMEGA_VALUE);
  u8_fe k1 = u8_fe_new(K1_VALUE);
  u8_fe k2 = u8_fe_new(K2_VALUE);

  size_t n = c->num_constraints;

  // step 2: create sigmas
  u8_fe *sigma_1 = (u8_fe *)malloc(n * sizeof(u8_fe));
  u8_fe *sigma_2 = (u8_fe *)malloc(n * sizeof(u8_fe));
  u8_fe *sigma_3 = (u8_fe *)malloc(n * sizeof(u8_fe));
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_prove\n");
    exit(EXIT_FAILURE);
  }

  copy_constraints_to_roots(p, c->c_a, n, sigma_1);
  copy_constraints_to_roots(p, c->c_b, n, sigma_2);
  copy_constraints_to_roots(p, c->c_c, n, sigma_3);

  // step 3: create polynomials
  poly f_a_x = interpolate_at_h(p, a->a, n);
  poly f_b_x = interpolate_at_h(p, a->b, n);
  poly f_c_x = interpolate_at_h(p, a->c, n);

  poly q_o_x = interpolate_at_h(p, c->q_o, n);
  poly q_m_x = interpolate_at_h(p, c->q_m, n);
  poly q_l_x = interpolate_at_h(p, c->q_l, n);
  poly q_r_x = interpolate_at_h(p, c->q_r, n);
  poly q_c_x = interpolate_at_h(p, c->q_c, n);

  poly s_sigma_1 = interpolate_at_h(p, sigma_1, n);
  poly s_sigma_2 = interpolate_at_h(p, sigma_2, n);
  poly s_sigma_3 = interpolate_at_h(p, sigma_3, n);

  // free sigma array as they are no longer needed
  free(sigma_1);
  free(sigma_2);
  free(sigma_3);

  // step 4: round 1 - evaluate a(x), b(x), c(x) at s
  u8_fe b1 = rand[0], b2 = rand[1], b3 = rand[2], b4 = rand[3], b5 = rand[4], b6 = rand[5];

  // a_x = (b2 + b1*x) * z_h_x + f_a_x
  u8_fe a_blinding_coeffs[] = {b2, b1};
  poly a_blinding_poly  = poly_new(a_blinding_coeffs, 2);
  poly a_x_blinded = poly_mul(&a_blinding_poly, &p->z_h_x);
  poly a_x = poly_add(&a_x_blinded, &f_a_x);

  // b_x = (b4 + b3*x) * z_h_x + f_b_x
  u8_fe b_blinding_coeffs[] = {b4, b3};
  poly b_blinding_poly  = poly_new(b_blinding_coeffs, 2);
  poly b_x_blinded = poly_mul(&b_blinding_poly, &p->z_h_x);
  poly b_x = poly_add(&b_x_blinded, &f_b_x);

  // c_x = (b6 + b5*x) * z_h_x + f_c_x
  u8_fe c_blinding_coeffs[] = {b6, b5};
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
  u8_fe b7 = rand[6], b8 = rand[7], b9 = rand[8];

  // initialize accumulator vector with 1
  u8_fe *acc = (u8_fe *)malloc((n + 1) *sizeof(u8_fe));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    exit(EXIT_FAILURE);
  }
  acc[0] = u8_fe_new(1);

  for (size_t i = 1; i <= n; i++) {
    // assignments at positoin i - 1
    u8_fe as_a = a->a[i - 1];
    u8_fe as_b = a->b[i - 1];
    u8_fe as_c = a->c[i - 1];

    // omega_pow = omega^(i-1)
    u8_fe omega_pow = u8_fe_pow(omega, i - 1);

    // denominator and numerator computations
    u8_fe denom = u8_fe_mul(
        u8_fe_mul(
            u8_fe_add(as_a, u8_fe_add(u8_fe_mul(beta, omega_pow), gamma)),
            u8_fe_add(as_b, u8_fe_add(u8_fe_mul(beta, u8_fe_mul(k1, omega_pow)), gamma))
                  ),
        u8_fe_add(as_c, u8_fe_add(u8_fe_mul(beta, u8_fe_mul(k2, omega_pow)), gamma))
			    );

    // s_sigma evaluations at omega_pow
    u8_fe s_sigma_1_eval = poly_eval(&s_sigma_1, omega_pow);
    u8_fe s_sigma_2_eval = poly_eval(&s_sigma_2, omega_pow);
    u8_fe s_sigma_3_eval = poly_eval(&s_sigma_3, omega_pow);

    u8_fe numer = u8_fe_mul(
        u8_fe_mul(
            u8_fe_add(as_a, u8_fe_add(u8_fe_mul(beta, s_sigma_1_eval), gamma)),
            u8_fe_add(as_b, u8_fe_add(u8_fe_mul(beta, s_sigma_2_eval), gamma))
                  ),
        u8_fe_add(as_c, u8_fe_add(u8_fe_mul(beta, s_sigma_3_eval), gamma))
                            );

    u8_fe fraction = u8_fe_div(denom, numer);
    acc[i] = u8_fe_mul(acc[i - 1], fraction);
  }

  // interpolate accumualtor vector
  poly acc_x = interpolate_at_h(p, acc, n + 1);
  free(acc);

  // check that acc_x evaluated at omega^n is 1
  u8_fe omega_n = u8_fe_pow(omega, n);
  u8_fe acc_x_eval = poly_eval(&acc_x, omega_n);
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

  // Step 10: Compute opening proofs w_z_x and w_z_omega_x
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

#endif // PLONK_H
