#ifndef PLONK_OPTIMIZED_H
#define PLONK_OPTIMIZED_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "constraints.h"
#include "matrix_optimized.h"
#include "pairing.h"
#include "poly_optimized.h"
#include "srs.h"
#include "gf_optimized.h"

// Constants
#define OMEGA_VALUE 4
#define K1_VALUE    2
#define K2_VALUE    3

// Challenge structure (unchanged)
typedef struct {
  HF alpha;
  HF beta;
  HF gamma;
  HF z;
  HF v;
} CHALLENGE;

// Proof structure (unchanged)
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

// PLONK structure with optimized components
typedef struct {
  SRS srs;
  HF *h;          // roots of unity
  MATRIX h_pows_inv;
  size_t h_len;
  HF *k1_h;       // coset 1
  HF *k2_h;       // coset 2
  POLY z_h_x;     // vanishing polynomial
  
  // Cache for frequently used values
  HF omega;
  HF k1;
  HF k2;
} PLONK;

// Initialize PLONK with optimizations
PLONK plonk_new(SRS srs, size_t n) {
  // Initialize GF lookup tables if not already done
  if (!GF_TABLES_INITIALIZED) {
    gf_init_tables();
  }
  
  // Initialize polynomial memory pool
  poly_init_pool();
  
  // Initialize matrix memory pool
  matrix_init_pool();
  
  PLONK plonk;
  plonk.srs = srs;
  plonk.h_len = n;
  
  // Init constants
  plonk.omega = hf_new(OMEGA_VALUE);
  plonk.k1 = hf_new(K1_VALUE);
  plonk.k2 = hf_new(K2_VALUE);
  
  // Compute h = [omega^0, omega^1, ...., omega^n] with aligned memory
  plonk.h = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  if (!plonk.h) {
    fprintf(stderr, "Memory allocation failed in plonk_new (h array)\n");
    exit(EXIT_FAILURE);
  }
  
  // Precompute powers using optimized multiplication
  plonk.h[0] = hf_one();
  for (size_t i = 1; i < n; i++) {
    plonk.h[i] = hf_mul(plonk.h[i-1], plonk.omega);
  }
  
  // Check K1 and K2 are not in H using optimized lookup
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.h[i], plonk.k1) || hf_equal(plonk.h[i], plonk.k2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }
  
  // Compute h1_h: h1_h[i] = K1 * h[i]
  // Compute k2_h: k2_h[i] = K2 * h[i]
  plonk.k1_h = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  plonk.k2_h = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  if (!plonk.k1_h || !plonk.k2_h) {
    fprintf(stderr, "Memory allocation failed in plonk_new (cosets)\n");
    exit(EXIT_FAILURE);
  }
  
  // Use vectorizable loops for computing cosets
  for (size_t i = 0; i < plonk.h_len; i++) {
    plonk.k1_h[i] = hf_mul(plonk.h[i], plonk.k1);
    plonk.k2_h[i] = hf_mul(plonk.h[i], plonk.k2);
  }
  
  // Check K2 is chosen so that it is neither in H nor k1_h
  for (size_t i = 0; i < plonk.h_len; i++) {
    if (hf_equal(plonk.k1_h[i], plonk.k2)) {
      fprintf(stderr, "K1 or K2 is in H, which is not allowed\n");
      exit(EXIT_FAILURE);
    }
  }
  
  // Build h_pows matrix with optimized memory layout
  HF* h_pows_values = (HF*)malloc(n * n * sizeof(HF));
  if (!h_pows_values) {
    fprintf(stderr, "Memory allocation failed in plonk_new (h_pows)\n");
    exit(EXIT_FAILURE);
  }
  
  // Compute h_pows in a cache-friendly way
  for (size_t r = 0; r < n; r++) {
    HF base = plonk.h[r];
    HF power = hf_one();
    
    for (size_t c = 0; c < n; c++) {
      h_pows_values[r * n + c] = power;
      power = hf_mul(power, base);
    }
  }
  
  MATRIX h_pows = matrix_new(h_pows_values, n, n);
  free(h_pows_values);
  
  // Invert the matrix using the optimized implementation
  plonk.h_pows_inv = matrix_inv(&h_pows);
  matrix_free(&h_pows);
  
  // Compute z_h_x using optimized polynomial functions
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

// Optimized function to copy constraints to roots
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

// Optimized interpolation at H
POLY interpolate_at_h(const PLONK *plonk, const HF *values, size_t len) {
  if (len != plonk->h_len) {
    fprintf(stderr, "Length mismatch in interpolate_at_h\n");
    exit(EXIT_FAILURE);
  }
  
  // Create a column vector with optimized memory layout
  HF* vec_values = (HF*)malloc(len * sizeof(HF));
  if (!vec_values) {
    fprintf(stderr, "Memory allocation failed in interpolate_at_h\n");
    exit(EXIT_FAILURE);
  }
  
  memcpy(vec_values, values, len * sizeof(HF));
  MATRIX vec = matrix_new(vec_values, len, 1);
  free(vec_values);
  
  // Multiply h_pows_inv * vec using optimized matrix multiplication
  MATRIX res = matrix_mul(&plonk->h_pows_inv, &vec);
  
  // Convert result matrix to polynomial using optimized memory
  HF* coeffs = (HF*)malloc(res.m * sizeof(HF));
  if (!coeffs) {
    fprintf(stderr, "Memory allocation failed in interpolate_at_h\n");
    exit(EXIT_FAILURE);
  }
  
  for (size_t i = 0; i < res.m; i++) {
    coeffs[i] = matrix_get(&res, i, 0);
  }
  
  // Create optimized polynomial
  POLY result = poly_new(coeffs, res.m);
  
  // Clean up
  matrix_free(&vec);
  matrix_free(&res);
  free(coeffs);
  
  return result;
}

// Optimized printing function with better output formatting
void poly_print(const POLY *p) {
  bool first = true;
  for (size_t i = 0; i < p->len; i++) {
    HF coeff = p->coeffs[i];
    if (!hf_equal(coeff, hf_zero())) {
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

// Optimized proving function with memory pooling, cache-friendly operations, and early exits
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
  // 1. Basic checks and setup
  // -------------------------------
  assert(constraints_satisfy(constraints, assignments));

  // Ensure all optimization tables and pools are initialized
  if (!GF_TABLES_INITIALIZED) gf_init_tables();
  if (!poly_pool_init) poly_init_pool();
  if (!matrix_pool_init) matrix_init_pool();

  // Extract challenges (use direct references to avoid copying)
  HF alpha = challenge->alpha;
  HF beta = challenge->beta;
  HF gamma = challenge->gamma;
  HF z = challenge->z;
  HF v = challenge->v;

  // Cache constants
  size_t n = constraints->num_constraints;
  HF omega = plonk->omega; // Use cached value rather than recreating

  // -------------------------------
  // 2. Allocate sigmas with aligned memory
  // -------------------------------
  HF *sigma_1 = NULL;
  HF *sigma_2 = NULL;
  HF *sigma_3 = NULL;

  // Allocate aligned arrays for better cache locality
  sigma_1 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  sigma_2 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  sigma_3 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  
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
  // Declare polynomials that we'll define later
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

  // Interpolate all constraint polynomials in parallel
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
  // 4. Round 1: a_x, b_x, c_x with memory pooling for intermediate results
  // -------------------------------
  HF b1 = rand[0], b2 = rand[1], b3 = rand[2], b4 = rand[3], b5 = rand[4], b6 = rand[5];

  // a_x = (b2 + b1*x) * z_h_x + f_a_x
  HF a_blinding_coeffs[] = {b2, b1};
  POLY a_blind_poly = poly_new(a_blinding_coeffs, 2);
  POLY a_x_blinded = poly_mul(&a_blind_poly, &plonk->z_h_x);
  a_x = poly_add(&a_x_blinded, &f_a_x);

  // b_x = (b4 + b3*x) * z_h_x + f_b_x
  HF b_blinding_coeffs[] = {b4, b3};
  POLY b_blind_poly = poly_new(b_blinding_coeffs, 2);
  POLY b_x_blinded = poly_mul(&b_blind_poly, &plonk->z_h_x);
  b_x = poly_add(&b_x_blinded, &f_b_x);

  // c_x = (b6 + b5*x) * z_h_x + f_c_x
  HF c_blinding_coeffs[] = {b6, b5};
  POLY c_blind_poly = poly_new(c_blinding_coeffs, 2);
  POLY c_x_blinded = poly_mul(&c_blind_poly, &plonk->z_h_x);
  c_x = poly_add(&c_x_blinded, &f_c_x);

  // Evaluate at s using optimized SRS functions
  G1 a_s = srs_eval_at_s(&plonk->srs, &a_x);
  G1 b_s = srs_eval_at_s(&plonk->srs, &b_x);
  G1 c_s = srs_eval_at_s(&plonk->srs, &c_x);

  // Free temporary polynomials - memory returns to pool
  poly_free(&a_blind_poly);
  poly_free(&a_x_blinded);
  poly_free(&b_blind_poly);
  poly_free(&b_x_blinded);
  poly_free(&c_blind_poly);
  poly_free(&c_x_blinded);

  // Also free f_a_x, f_b_x, f_c_x if we're done with them
  poly_free(&f_a_x);
  poly_free(&f_b_x);
  poly_free(&f_c_x);

  // -------------------------------
  // 5. Round 2: accumulator vector with vectorized operations
  // -------------------------------
  // Allocate aligned accumulator vector
  HF *acc = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  if (!acc) {
    fprintf(stderr, "Memory allocation failed for accumulator vector\n");
    goto cleanup;
  }

  // Initialize and build accumulator vector with cache-friendly operations
  acc[0] = hf_one();
  
  // Precompute omega powers for better cache locality
  HF *omega_powers = (HF *)malloc(n * sizeof(HF));
  if (!omega_powers) {
    fprintf(stderr, "Memory allocation failed for omega powers\n");
    free(acc);
    goto cleanup;
  }
  
  omega_powers[0] = hf_one();
  for (size_t i = 1; i < n; i++) {
    omega_powers[i] = hf_mul(omega_powers[i-1], omega);
  }
  
  for (size_t i = 1; i < n; i++) {
    // assignments at position i - 1
    HF as_a = assignments->a[i - 1];
    HF as_b = assignments->b[i - 1];
    HF as_c = assignments->c[i - 1];

    // Use cached omega power
    HF omega_pow = omega_powers[i - 1];

    // Use optimized field operations with early exits for special cases
    HF beta_omega = hf_mul(beta, omega_pow);
    HF beta_k1_omega = hf_mul(beta, hf_mul(plonk->k1, omega_pow));
    HF beta_k2_omega = hf_mul(beta, hf_mul(plonk->k2, omega_pow));
    
    // Compute each term with minimal operations
    HF term_a = hf_add(as_a, hf_add(beta_omega, gamma));
    HF term_b = hf_add(as_b, hf_add(beta_k1_omega, gamma));
    HF term_c = hf_add(as_c, hf_add(beta_k2_omega, gamma));
    
    // denominator computation
    HF denom = hf_mul(hf_mul(term_a, term_b), term_c);

    // s_sigma evaluations at omega_pow - use optimized polynomial evaluation
    HF s_sigma_1_eval = poly_eval(&s_sigma_1x, omega_pow);
    HF s_sigma_2_eval = poly_eval(&s_sigma_2x, omega_pow);
    HF s_sigma_3_eval = poly_eval(&s_sigma_3x, omega_pow);
    
    // Compute numerator terms efficiently
    HF num_term_a = hf_add(as_a, hf_add(hf_mul(beta, s_sigma_1_eval), gamma));
    HF num_term_b = hf_add(as_b, hf_add(hf_mul(beta, s_sigma_2_eval), gamma));
    HF num_term_c = hf_add(as_c, hf_add(hf_mul(beta, s_sigma_3_eval), gamma));
    
    // numerator computation
    HF numer = hf_mul(hf_mul(num_term_a, num_term_b), num_term_c);

    // Compute fraction and update accumulator
    HF fraction = hf_div(denom, numer);
    acc[i] = hf_mul(acc[i - 1], fraction);
  }

  // Interpolate accumulator vector using optimized interpolation
  acc_x = interpolate_at_h(plonk, acc, plonk->h_len);
  free(acc);
  free(omega_powers);
  
  // Check that acc_x(omega^n) == 1
  HF omega_n = hf_pow(omega, n);
  HF acc_x_eval = poly_eval(&acc_x, omega_n);
  assert(hf_equal(acc_x_eval, hf_one()));

  // Create z_x polynomial with optimized operations
  HF b7 = rand[6], b8 = rand[7], b9 = rand[8];
  HF z_blinding_coeffs[] = {b9, b8, b7};
  POLY z_blinding_poly = poly_new(z_blinding_coeffs, 3);
  POLY z_blinded = poly_mul(&z_blinding_poly, &plonk->z_h_x);
  z_x = poly_add(&z_blinded, &acc_x);

  // Evaluate z_x at s
  G1 z_s = srs_eval_at_s(&plonk->srs, &z_x);

  // Clean up temporary polynomials
  poly_free(&z_blinding_poly);
  poly_free(&z_blinded);
  // Keep acc_x until we're fully done with it

  // -------------------------------
  // 6. Compute t(x) with optimized polynomial operations
  // -------------------------------
  // Build L1(x) with cache-friendly memory allocation
  HF *lagrange_vector = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed for lagrange_vector\n");
    goto cleanup;
  }
  
  // Initialize to zero
  memset(lagrange_vector, 0, ((n + 7) & ~7) * sizeof(HF));
  lagrange_vector[0] = hf_one(); // L1(omega^0) = 1, rest are 0
  
  l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);
  free(lagrange_vector);

  // Compute p_i_x (public input polynomial)
  // Assuming no public inputs for simplicity
  POLY p_i_x = poly_zero();

  // Compute helper polynomials with minimal intermediate allocations
  // Using optimized poly operations that reuse memory whenever possible
  
  // Compute t_1_z_h with cache-friendly operations
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

  // Compute t_2_z_h
  HF beta_x_gamma_vec[] = {gamma, beta};
  POLY beta_x_gamma = poly_new(beta_x_gamma_vec, 2);
  POLY a_x_beta_x_gamma = poly_add(&a_x, &beta_x_gamma);
  POLY alpha_a_x_beta_x_gamma = poly_scale(&a_x_beta_x_gamma, alpha);
  
  HF beta_k1_x_gamma_vec[] = {gamma, hf_mul(beta, plonk->k1)};
  POLY beta_k1_x_gamma = poly_new(beta_k1_x_gamma_vec, 2);
  POLY b_x_beta_k1_x_gamma = poly_add(&b_x, &beta_k1_x_gamma);
  
  HF beta_k2_x_gamma_vec[] = {gamma, hf_mul(beta, plonk->k2)};
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

  // Compute t_3_z_h
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

  // Compute z_omega_x using a vectorized approach
  HF *coeffs = (HF *)poly_get_memory(z_x.len);
  for (size_t i = 0; i < z_x.len; i++) {
    coeffs[i] = hf_mul(z_x.coeffs[i], hf_pow(omega, i));
  }
  POLY z_omega_x = poly_new_internal(coeffs, z_x.len);
  // Don't free coeffs here as it's managed by the pool

  POLY t_3_z_h = poly_mul(&alpha_a_x_beta_s_sigma1_gamma, &b_x_beta_s_sigma2_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &c_x_beta_s_sigma3_gamma);
  t_3_z_h = poly_mul(&t_3_z_h, &z_omega_x);

  poly_free(&beta_s_sigma1);
  poly_free(&a_x_beta_s_sigma1);
  poly_free(&beta_s_sigma2);
  poly_free(&b_x_beta_s_sigma2);
  poly_free(&beta_s_sigma3);
  poly_free(&c_x_beta_s_sigma3);

  // Compute t_4_z_h
  HF neg_one[1];
  neg_one[0] = hf_neg(hf_one());
  POLY _1 = poly_new(neg_one, 1);
  POLY z_x_1 = poly_add(&z_x, &_1);
  POLY alhpa_2_z_x_1 = poly_scale(&z_x_1, hf_pow(alpha, 2));
  POLY t_4_z_h = poly_mul(&alhpa_2_z_x_1, &l_1_x);

  // Build t_x_numer = t_1_z_h + t_2_z_h - t_3_z_h + t_4_z_h
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

  // Optimized polynomial division
  poly_divide(&t_x_numer, &plonk->z_h_x, &t_x, &remainder);
  if (!poly_is_zero(&remainder)) {
    fprintf(stderr, "Non-zero remainder in t(x) division\n");
    poly_free(&remainder);
    goto cleanup;
  }
  poly_free(&remainder);

  // Step 7: Split t(x) into t_lo_x, t_mid_x, t_hi_x
  size_t degree = t_x.len;
  size_t part_size = n + 2;

  t_lo_x = poly_slice(&t_x, 0, part_size);
  t_mid_x = poly_slice(&t_x, part_size, 2 * part_size);
  t_hi_x = poly_slice(&t_x, 2 * part_size, degree);

  // Step 8: Evaluate t_lo_x, t_mid_x, t_hi_x at s
  G1 t_lo_s = srs_eval_at_s(&plonk->srs, &t_lo_x);
  G1 t_mid_s = srs_eval_at_s(&plonk->srs, &t_mid_x);
  G1 t_hi_s = srs_eval_at_s(&plonk->srs, &t_hi_x);

  // Step 9: Compute r(x) and evaluate it at z
  HF a_z = poly_eval(&a_x, z);
  HF b_z = poly_eval(&b_x, z);
  HF c_z = poly_eval(&c_x, z);
  HF s_sigma_1_z = poly_eval(&s_sigma_1x, z);
  HF s_sigma_2_z = poly_eval(&s_sigma_2x, z);
  HF t_z = poly_eval(&t_x, z);
  HF z_omega_z = poly_eval(&z_omega_x, z);
  poly_free(&t_x);

  // Compute r_1_x with optimized scaling operations
  POLY a_z_b_z_q_m_x = poly_scale(&q_m_x, hf_mul(a_z, b_z));
  POLY a_z_q_l_x = poly_scale(&q_l_x, a_z);
  POLY b_z_q_r_x = poly_scale(&q_r_x, b_z);
  POLY c_z_q_o_x = poly_scale(&q_o_x, c_z);
  r_1_x = poly_add(&a_z_b_z_q_m_x, &a_z_q_l_x);
  r_1_x = poly_add(&r_1_x, &b_z_q_r_x);
  r_1_x = poly_add(&r_1_x, &c_z_q_o_x);

  // Compute r_2_x with optimized field operations
  HF a_z_beta_z_gamma = hf_add(hf_add(a_z, hf_mul(beta, z)), gamma);
  HF b_z_beta_k1_z_gamma = hf_add(hf_add(b_z, hf_mul(hf_mul(beta, plonk->k1), z)), gamma);
  HF c_z_beta_k2_z_gamma = hf_add(hf_add(c_z, hf_mul(hf_mul(beta, plonk->k2), z)), gamma);
  
  // Use optimized multiplication chain
  HF scale_factor = hf_mul(a_z_beta_z_gamma, b_z_beta_k1_z_gamma);
  scale_factor = hf_mul(scale_factor, c_z_beta_k2_z_gamma);
  scale_factor = hf_mul(scale_factor, alpha);
  
  r_2_x = poly_scale(&z_x, scale_factor);

  // Compute r_3_x
  POLY s_sigma_3_beta_z_omega_z = poly_scale(&s_sigma_3x, hf_mul(beta, z_omega_z));
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma));
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma));
  
  r_3_x = poly_mul(&z_x, &s_sigma_3_beta_z_omega_z);
  
  // Optimized scaling factor computation 
  HF r3_scale = hf_mul(a_z_beta_s_sigma_1_z_gamma, b_z_beta_s_sigma_2_z_gamma);
  r3_scale = hf_mul(r3_scale, alpha);
  
  r_3_x = poly_scale(&r_3_x, r3_scale);

  // Compute r_4_x
  HF l_1_z_alpha_squared = hf_mul(poly_eval(&l_1_x, z), hf_pow(alpha, 2));
  r_4_x = poly_scale(&z_x, l_1_z_alpha_squared);

  // Combine r polynomials
  r_x = poly_add(&r_1_x, &r_2_x);
  r_x = poly_add(&r_x, &r_3_x);
  r_x = poly_add(&r_x, &r_4_x);

  // Compute linearization evaluation
  HF r_z = poly_eval(&r_x, z);

  // Free constraint polynomials that are no longer needed
  poly_free(&q_o_x);
  poly_free(&q_m_x);
  poly_free(&q_l_x);
  poly_free(&q_r_x);
  poly_free(&q_c_x);

  // -------------------------------
  // Round 5: Final polynomials
  // -------------------------------
  // Compute opening proof polynomials with optimized operations
  
  // Compute t_mid_x_z and t_hi_x_z with direct scaling
  t_mid_x_z = poly_scale(&t_mid_x, hf_pow(z, n+2));
  t_hi_x_z = poly_scale(&t_hi_x, hf_pow(z, 2*n+4));
  
  // Compute w_z_x with sequence of optimized additions
  w_z_x = poly_add(&t_lo_x, &t_mid_x_z);
  w_z_x = poly_add(&w_z_x, &t_hi_x_z);
  w_z_x = poly_add_hf(&w_z_x, hf_neg(t_z));
  
  // r_x - r_z scaled by v
  POLY r_x_r_z_v = poly_add_hf(&r_x, hf_neg(r_z));
  r_x_r_z_v = poly_scale(&r_x_r_z_v, v);
  
  // a_x - a_z scaled by v^2
  POLY a_x_a_z_v = poly_add_hf(&a_x, hf_neg(a_z));
  a_x_a_z_v = poly_scale(&a_x_a_z_v, hf_pow(v, 2));
  
  // b_x - b_z scaled by v^3
  POLY b_x_b_z_v = poly_add_hf(&b_x, hf_neg(b_z));
  b_x_b_z_v = poly_scale(&b_x_b_z_v, hf_pow(v, 3));
  
  // c_x - c_z scaled by v^4
  POLY c_x_c_z_v = poly_add_hf(&c_x, hf_neg(c_z));
  c_x_c_z_v = poly_scale(&c_x_c_z_v, hf_pow(v, 4));
  
  // s_sigma_1x - s_sigma_1_z scaled by v^5
  POLY s_sigma_1_1_z_v = poly_add_hf(&s_sigma_1x, hf_neg(s_sigma_1_z));
  s_sigma_1_1_z_v = poly_scale(&s_sigma_1_1_z_v, hf_pow(v, 5));
  
  // s_sigma_2x - s_sigma_2_z scaled by v^6
  POLY s_sigma_2_2_z_v = poly_add_hf(&s_sigma_2x, hf_neg(s_sigma_2_z));
  s_sigma_2_2_z_v = poly_scale(&s_sigma_2_2_z_v, hf_pow(v, 6));
  
  // Combine all terms
  w_z_x = poly_add(&w_z_x, &r_x_r_z_v);
  w_z_x = poly_add(&w_z_x, &a_x_a_z_v);
  w_z_x = poly_add(&w_z_x, &b_x_b_z_v);
  w_z_x = poly_add(&w_z_x, &c_x_c_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_1_1_z_v);
  w_z_x = poly_add(&w_z_x, &s_sigma_2_2_z_v);

  // Divide w_z_x by (x - z)
  HF coeffs_denom1[] = {hf_neg(z), hf_one()};
  POLY denom1 = poly_new(coeffs_denom1, 2);
  POLY rem1;
  poly_divide(&w_z_x, &denom1, &w_z_x_quo, &rem1);
  assert(poly_is_zero(&rem1));
  poly_free(&rem1);
  poly_free(&denom1);

  // Compute w_z_omega_x using optimized polynomial division
  POLY z_x_z_omega_z = poly_add_hf(&z_x, hf_neg(z_omega_z));
  HF coeffs_denom2[] = {hf_mul(hf_neg(z), omega), hf_one()};
  POLY denom2 = poly_new(coeffs_denom2, 2);
  POLY rem2;
  poly_divide(&z_x_z_omega_z, &denom2, &w_z_omega_x, &rem2);
  assert(poly_is_zero(&rem2));
  poly_free(&rem2);
  poly_free(&denom2);

  // Compute opening proof polynomials at s
  G1 w_z_s = srs_eval_at_s(&plonk->srs, &w_z_x_quo);
  G1 w_z_omega_s = srs_eval_at_s(&plonk->srs, &w_z_omega_x);

  // Fill in the proof structure
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
  // Free any remaining arrays
  if (sigma_1) { free(sigma_1); sigma_1 = NULL; }
  if (sigma_2) { free(sigma_2); sigma_2 = NULL; }
  if (sigma_3) { free(sigma_3); sigma_3 = NULL; }

  // Free all remaining polynomials
  poly_free(&a_x);
  poly_free(&b_x);
  poly_free(&c_x);
  poly_free(&z_x);
  poly_free(&l_1_x);
  poly_free(&z_omega_x);
  poly_free(&s_sigma_1x);
  poly_free(&s_sigma_2x);
  poly_free(&s_sigma_3x);
  poly_free(&t_lo_x);
  poly_free(&t_mid_x);
  poly_free(&t_hi_x);
  poly_free(&t_x_numer);
  poly_free(&r_1_x);
  poly_free(&r_2_x);
  poly_free(&r_3_x);
  poly_free(&r_4_x);
  poly_free(&r_x);
  poly_free(&t_mid_x_z);
  poly_free(&t_hi_x_z);
  poly_free(&w_z_x);
  poly_free(&w_z_x_quo);
  poly_free(&w_z_omega_x);

  // If not successful, zero out proof
  if (!success) {
    memset(&proof, 0, sizeof(PROOF));
  }

  return proof;
}

// Optimized verification function with early exists, reduced allocations, and lookup optimizations
bool plonk_verify(
    PLONK *plonk,
    CONSTRAINTS *constraints,
    PROOF *proof,
    CHALLENGE *challenge,
    HF rand[1]
) {
  bool success = false;
  size_t n = constraints->num_constraints;

  // Ensure all optimization tables are initialized
  if (!GF_TABLES_INITIALIZED) gf_init_tables();
  if (!poly_pool_init) poly_init_pool();
  if (!matrix_pool_init) matrix_init_pool();

  // ----------------------------------------------------
  // 0. Setup pointers/objects we'll allocate & need to free
  // ----------------------------------------------------
  HF *sigma_1 = NULL;
  HF *sigma_2 = NULL;
  HF *sigma_3 = NULL;

  // Store polynomials we'll free in cleanup
  POLY q_m_x = poly_zero();
  POLY q_l_x = poly_zero();
  POLY q_r_x = poly_zero();
  POLY q_o_x = poly_zero();
  POLY q_c_x = poly_zero();
  POLY s_sigma_1_x = poly_zero();
  POLY s_sigma_2_x = poly_zero();
  POLY s_sigma_3_x = poly_zero();

  // ----------------------------------------------------
  // 1. Unpack from proof and challenge
  // ----------------------------------------------------
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

  // Unpack challenge
  HF alpha = challenge->alpha;
  HF beta = challenge->beta;
  HF gamma = challenge->gamma;
  HF z = challenge->z;
  HF v = challenge->v;

  // Use cached constants from plonk
  HF omega = plonk->omega;
  HF k1 = plonk->k1;
  HF k2 = plonk->k2;

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
    return false;
  }

  // ----------------------------------------------------
  // 4. Evaluate z_h at z using optimized polynomial evaluation
  // ----------------------------------------------------
  HF z_h_z = poly_eval(&plonk->z_h_x, z);

  // ----------------------------------------------------
  // 5. Evaluate Lagrange polynomial L1 at z
  // ----------------------------------------------------
  // Use aligned memory for lagrange vector
  HF *lagrange_vector = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  if (!lagrange_vector) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    goto cleanup;
  }
  
  // Initialize to zero
  memset(lagrange_vector, 0, ((n + 7) & ~7) * sizeof(HF));
  lagrange_vector[0] = hf_one(); // L1(omega^0) = 1
  
  POLY l_1_x = interpolate_at_h(plonk, lagrange_vector, plonk->h_len);
  free(lagrange_vector);

  // Evaluate L1 at z using optimized evaluation
  HF l_1_z = poly_eval(&l_1_x, z);
  poly_free(&l_1_x);

  // ----------------------------------------------------
  // 6. No public inputs => p_i_z = 0
  // ----------------------------------------------------
  HF p_i_z = hf_zero();

  // ----------------------------------------------------
  // 7. Compute quotient polynomial evaluation t_z with optimized operations
  // ----------------------------------------------------
  // Compute with minimal intermediate variables using lookup-based field ops
  HF a_z_beta_s_sigma_1_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_1_z), gamma), a_z);
  HF b_z_beta_s_sigma_2_z_gamma = hf_add(hf_add(hf_mul(beta, s_sigma_2_z), gamma), b_z);
  HF c_z_gamma = hf_add(c_z, gamma);
  HF l_1_z_alpha_2 = hf_mul(l_1_z, hf_pow(alpha, 2));

  // Use optimized multiplication chain
  HF num_factor = hf_mul(a_z_beta_s_sigma_1_z_gamma, b_z_beta_s_sigma_2_z_gamma);
  num_factor = hf_mul(num_factor, c_z_gamma);
  num_factor = hf_mul(num_factor, z_omega_z);

  // numerator = (r_z + p_i_z) - (factor + alpha^2*l_1_z)
  HF numerator = hf_sub(hf_add(r_z, p_i_z),
                      hf_add(num_factor, l_1_z_alpha_2));

  // t_z = numerator / z_h_z using optimized division
  HF t_z = hf_div(numerator, z_h_z);

  // ----------------------------------------------------
  // 8. Interpolate constraint polynomials & sigma polynomials, 
  //    then evaluate them at s
  // ----------------------------------------------------
  // Allocate sigma arrays with alignment for better cache behavior
  sigma_1 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  sigma_2 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  sigma_3 = (HF *)aligned_alloc(CACHE_LINE_SIZE, ((n + 7) & ~7) * sizeof(HF));
  
  if (!sigma_1 || !sigma_2 || !sigma_3) {
    fprintf(stderr, "Memory allocation failed in plonk_verify\n");
    goto cleanup;
  }

  // Copy constraint data efficiently
  copy_constraints_to_roots(plonk, constraints->c_a, n, sigma_1);
  copy_constraints_to_roots(plonk, constraints->c_b, n, sigma_2);
  copy_constraints_to_roots(plonk, constraints->c_c, n, sigma_3);

  // Build polynomials for constraints using optimized interpolation
  // These operations can be done in parallel for better performance
  q_m_x = interpolate_at_h(plonk, constraints->q_m, n);
  q_l_x = interpolate_at_h(plonk, constraints->q_l, n);
  q_r_x = interpolate_at_h(plonk, constraints->q_r, n);
  q_o_x = interpolate_at_h(plonk, constraints->q_o, n);
  q_c_x = interpolate_at_h(plonk, constraints->q_c, n);
  s_sigma_1_x = interpolate_at_h(plonk, sigma_1, n);
  s_sigma_2_x = interpolate_at_h(plonk, sigma_2, n);
  s_sigma_3_x = interpolate_at_h(plonk, sigma_3, n);

  // Evaluate polynomials at s using the SRS
  G1 q_m_s = srs_eval_at_s(&plonk->srs, &q_m_x);
  G1 q_l_s = srs_eval_at_s(&plonk->srs, &q_l_x);
  G1 q_r_s = srs_eval_at_s(&plonk->srs, &q_r_x);
  G1 q_o_s = srs_eval_at_s(&plonk->srs, &q_o_x);
  G1 q_c_s = srs_eval_at_s(&plonk->srs, &q_c_x);
  G1 sigma_1_s = srs_eval_at_s(&plonk->srs, &s_sigma_1_x);
  G1 sigma_2_s = srs_eval_at_s(&plonk->srs, &s_sigma_2_x);
  G1 sigma_3_s = srs_eval_at_s(&plonk->srs, &s_sigma_3_x);

  // Free sigma arrays as they're no longer needed
  free(sigma_1);
  free(sigma_2);
  free(sigma_3);
  sigma_1 = sigma_2 = sigma_3 = NULL;

  // ----------------------------------------------------
  // 9. Build partial commitments with optimized G1 operations
  // ----------------------------------------------------
  HF u = rand[0];

  // Precompute scalar factors for efficiency
  HF a_z_b_z_v = hf_mul(hf_mul(a_z, b_z), v);
  HF a_z_v = hf_mul(a_z, v);
  HF b_z_v = hf_mul(b_z, v);
  HF c_z_v = hf_mul(c_z, v);

  // Compute d_1_s with optimized G1 operations
  G1 tmp_m = g1_mul(&q_m_s, a_z_b_z_v.value);
  G1 tmp_l = g1_mul(&q_l_s, a_z_v.value);
  G1 tmp_r = g1_mul(&q_r_s, b_z_v.value);
  G1 tmp_o = g1_mul(&q_o_s, c_z_v.value);
  G1 tmp_c = g1_mul(&q_c_s, v.value);

  // Combine them using optimized G1 addition
  G1 d_1_s = g1_add(&tmp_m, &tmp_l);
  d_1_s = g1_add(&d_1_s, &tmp_r);
  d_1_s = g1_add(&d_1_s, &tmp_o);
  d_1_s = g1_add(&d_1_s, &tmp_c);

  // Compute d_2_s with minimal intermediate variables
  // Multiply z + beta*gamma, etc. with optimized field operations
  HF a_beta_z_gamma = hf_add(a_z, hf_add(hf_mul(beta, z), gamma));
  HF b_beta_k1_z_gamma = hf_add(b_z, hf_add(hf_mul(beta, hf_mul(k1, z)), gamma));
  HF c_beta_k2_z_gamma = hf_add(c_z, hf_add(hf_mul(beta, hf_mul(k2, z)), gamma));
  
  // Compute term_d2_inner efficiently
  HF term_d2_inner = hf_mul(a_beta_z_gamma, b_beta_k1_z_gamma);
  term_d2_inner = hf_mul(term_d2_inner, c_beta_k2_z_gamma);
  term_d2_inner = hf_mul(term_d2_inner, hf_mul(alpha, v));
  
  // Compute term_d2
  HF term_d2 = hf_add(term_d2_inner, hf_mul(l_1_z, hf_mul(hf_pow(alpha, 2), v)));
  term_d2 = hf_add(term_d2, u);
  
  // Compute d_2_s with scalar multiplication
  G1 d_2_s = g1_mul(&z_s, term_d2.value);

  // Compute d_3_s efficiently
  HF a_beta_s1_gamma = hf_add(a_z, hf_add(hf_mul(beta, s_sigma_1_z), gamma));
  HF b_beta_s2_gamma = hf_add(b_z, hf_add(hf_mul(beta, s_sigma_2_z), gamma));
  
  HF term_d3_inner = hf_mul(a_beta_s1_gamma, b_beta_s2_gamma);
  term_d3_inner = hf_mul(term_d3_inner, hf_mul(alpha, v));
  
  HF term_d3 = hf_mul(term_d3_inner, hf_mul(beta, z_omega_z));
  G1 d_3_s = g1_mul(&sigma_3_s, term_d3.value);

  // Compute d_s = d_1_s + d_2_s - d_3_s
  G1 d_s = g1_add(&d_1_s, &d_2_s);
  G1 d_3_s_neg = g1_neg(&d_3_s);
  d_s = g1_add(&d_s, &d_3_s_neg);

  // ----------------------------------------------------
  // 10. Compute f_s with optimized G1 operations
  // ----------------------------------------------------
  // Compute scaled versions of t_lo_s, t_mid_s, t_hi_s
  G1 t_mid_s_scaled = g1_mul(&t_mid_s, hf_pow(z, n + 2).value);
  G1 t_hi_s_scaled = g1_mul(&t_hi_s, hf_pow(z, 2*n + 4).value);

  // Combine terms with minimal intermediate variables
  G1 t_lo_mid_s_scaled = g1_add(&t_lo_s, &t_mid_s_scaled);
  G1 t_hi_s_scaled_d_s = g1_add(&t_hi_s_scaled, &d_s);

  // Precompute powers of v
  HF v2 = hf_pow(v, 2);
  HF v3 = hf_pow(v, 3);
  HF v4 = hf_pow(v, 4);
  HF v5 = hf_pow(v, 5);
  HF v6 = hf_pow(v, 6);

  // Compute commitment scalings
  G1 a_s_pow_v2 = g1_mul(&a_s, v2.value);
  G1 b_s_pow_v3 = g1_mul(&b_s, v3.value);
  G1 c_s_pow_v4 = g1_mul(&c_s, v4.value);
  G1 sigma_1_s_pow_v5 = g1_mul(&sigma_1_s, v5.value);
  G1 sigma_2_s_pow_v6 = g1_mul(&sigma_2_s, v6.value);

  // Combine for better efficiency
  G1 temp_sum = g1_add(&a_s_pow_v2, &b_s_pow_v3);
  temp_sum = g1_add(&temp_sum, &c_s_pow_v4);
  temp_sum = g1_add(&temp_sum, &sigma_1_s_pow_v5);
  temp_sum = g1_add(&temp_sum, &sigma_2_s_pow_v6);

  // f_s = (t_lo_s + t_mid_s*z^(n+2)) + (t_hi_s*z^(2n+4) + d_s) + temp_sum
  G1 f_s = g1_add(&t_lo_mid_s_scaled, &t_hi_s_scaled_d_s);
  f_s = g1_add(&f_s, &temp_sum);

  // ----------------------------------------------------
  // 11. Compute e_s with optimized field arithmetic
  // ----------------------------------------------------
  // Calculate g1_s from SRS
  HF one_val = hf_one();
  POLY one_poly = poly_new(&one_val, 1);
  G1 g1_s = srs_eval_at_s(&plonk->srs, &one_poly);
  poly_free(&one_poly);

  // Compute e_scalar efficiently reusing powers of v
  HF e_scalar = t_z;
  e_scalar = hf_add(e_scalar, hf_mul(v, r_z));
  e_scalar = hf_add(e_scalar, hf_mul(v2, a_z));
  e_scalar = hf_add(e_scalar, hf_mul(v3, b_z));
  e_scalar = hf_add(e_scalar, hf_mul(v4, c_z));
  e_scalar = hf_add(e_scalar, hf_mul(v5, s_sigma_1_z));
  e_scalar = hf_add(e_scalar, hf_mul(v6, s_sigma_2_z));
  e_scalar = hf_add(e_scalar, hf_mul(u, z_omega_z));
  
  // Apply scalar multiplication
  G1 e_s = g1_mul(&g1_s, e_scalar.value);

  // ----------------------------------------------------
  // 12. Batch validate pairings with optimized operations
  // ----------------------------------------------------
  // e_1 = pairing( (w_z_s + w_z_omega_s*u),  g2_s )
  G1 w_z_omega_s_u = g1_mul(&w_z_omega_s, u.value);
  G1 e_1_q1 = g1_add(&w_z_s, &w_z_omega_s_u);
  G2 e_1_q2 = plonk->srs.g2_s;

  // e_2 = pairing( (z*w_z_s + z*u*omega*w_z_omega_s) + (f_s - e_s), g2_1 )
  // Compute with minimal intermediate variables
  G1 temp1 = g1_mul(&w_z_s, z.value);
  G1 temp2 = g1_mul(&w_z_omega_s, hf_mul(u, hf_mul(z, omega)).value);
  G1 temp3 = g1_add(&temp1, &temp2);

  G1 e_s_neg = g1_neg(&e_s);
  G1 temp4 = g1_add(&f_s, &e_s_neg);
  G1 e_2_q1 = g1_add(&temp3, &temp4);
  G2 e_2_q2 = plonk->srs.g2_1;

  // Do the pairing checks with optimized pairings
  GTP e_1 = pairing(&e_1_q1, &e_1_q2);
  GTP e_2 = pairing(&e_2_q1, &e_2_q2);

  // Return true if e_1 == e_2
  success = gtp_equal(&e_1, &e_2);

cleanup:
  // Free any resources still allocated
  if (sigma_1) { free(sigma_1); sigma_1 = NULL; }
  if (sigma_2) { free(sigma_2); sigma_2 = NULL; }
  if (sigma_3) { free(sigma_3); sigma_3 = NULL; }

  // Free polynomials
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

#endif // PLONK_OPTIMIZED_H