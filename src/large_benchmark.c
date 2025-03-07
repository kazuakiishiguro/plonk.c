#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "plonk.h"

double get_time_ms() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec * 1000.0) + (tv.tv_usec / 1000.0);
}

void benchmark_large_plonk() {
  printf("Running large PLONK computation benchmark...\n");
  
  // Create the trusted setup - use the same parameters as in the test for now
  GF secret = f101(2); // the toxic waste
  size_t n = 6;       // Polynomial degree (must be large enough for h_len)
  size_t h_len = 4;   // Constraints
  
  double start, end;
  
  // Benchmark SRS creation
  start = get_time_ms();
  SRS srs = srs_create(secret, n);
  end = get_time_ms();
  printf("Large SRS creation: %.2f ms\n", end - start);
  
  // Benchmark PLONK initialization with large parameters
  start = get_time_ms();
  PLONK plonk = plonk_new(srs, h_len);
  end = get_time_ms();
  printf("Large PLONK initialization: %.2f ms\n", end - start);

  // Set up constraints for a larger circuit
  start = get_time_ms();
  CONSTRAINTS constraints;
  constraints.num_constraints = h_len;

  // Initialize q_m, q_l, q_r, q_o, q_c arrays
  constraints.q_m = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_l = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_r = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_o = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  constraints.q_c = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  
  if (!constraints.q_m || !constraints.q_l || !constraints.q_r || !constraints.q_o || !constraints.q_c) {
    fprintf(stderr, "Memory allocation failed in benchmark\n");
    exit(EXIT_FAILURE);
  }
  
  // Create a mix of multiplication and addition gates
  for (size_t i = 0; i < h_len; i++) {
    if (i % 3 == 0) {
      // Multiplication gate: a * b - c = 0
      constraints.q_m[i] = hf_one();
      constraints.q_l[i] = hf_zero();
      constraints.q_r[i] = hf_zero();
      constraints.q_o[i] = hf_neg(hf_one());
      constraints.q_c[i] = hf_zero();
    } else if (i % 3 == 1) {
      // Addition gate: a + b - c = 0
      constraints.q_m[i] = hf_zero();
      constraints.q_l[i] = hf_one();
      constraints.q_r[i] = hf_one();
      constraints.q_o[i] = hf_neg(hf_one());
      constraints.q_c[i] = hf_zero();
    } else {
      // Subtraction gate: a - b - c = 0
      constraints.q_m[i] = hf_zero();
      constraints.q_l[i] = hf_one();
      constraints.q_r[i] = hf_neg(hf_one());
      constraints.q_o[i] = hf_neg(hf_one());
      constraints.q_c[i] = hf_zero();
    }
  }

  // Copy constraints
  constraints.c_a = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_b = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));
  constraints.c_c = (COPY_OF *)malloc(constraints.num_constraints * sizeof(COPY_OF));

  if (!constraints.c_a || !constraints.c_b || !constraints.c_c) {
    fprintf(stderr, "Memory allocation failed in benchmark\n");
    exit(EXIT_FAILURE);
  }

  // Set up copy constraints
  for (size_t i = 0; i < h_len; i++) {
    constraints.c_a[i] = (COPY_OF){COPYOF_A, i + 1};
    constraints.c_b[i] = (COPY_OF){COPYOF_B, i + 1};
    constraints.c_c[i] = (COPY_OF){COPYOF_C, i + 1};
  }
  end = get_time_ms();
  printf("Large constraints setup: %.2f ms\n", end - start);

  // Define assignments that satisfy the constraints
  start = get_time_ms();
  ASSIGNMENTS assignments;
  assignments.len = constraints.num_constraints;
  assignments.a = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  assignments.b = (HF *)malloc(constraints.num_constraints * sizeof(HF));
  assignments.c = (HF *)malloc(constraints.num_constraints * sizeof(HF));

  if (!assignments.a || !assignments.b || !assignments.c) {
    fprintf(stderr, "Memory allocation failed in benchmark\n");
    exit(EXIT_FAILURE);
  }

  // Generate values for a and b, then compute c based on the gate type
  for (size_t i = 0; i < h_len; i++) {
    assignments.a[i] = hf_new(i + 3);  // Some arbitrary value
    assignments.b[i] = hf_new(i + 5);  // Some arbitrary value
    
    if (i % 3 == 0) {
      // Multiplication gate: c = a * b
      assignments.c[i] = hf_mul(assignments.a[i], assignments.b[i]);
    } else if (i % 3 == 1) {
      // Addition gate: c = a + b
      assignments.c[i] = hf_add(assignments.a[i], assignments.b[i]);
    } else {
      // Subtraction gate: c = a - b
      assignments.c[i] = hf_sub(assignments.a[i], assignments.b[i]);
    }
  }
  end = get_time_ms();
  printf("Large assignments setup: %.2f ms\n", end - start);

  // Random numbers for blinding
  HF rand[9] = {
    hf_new(7), hf_new(4), hf_new(11),
    hf_new(12), hf_new(16), hf_new(2),
    hf_new(14), hf_new(11), hf_new(7)
  };

  // Challenge values
  CHALLENGE challenge;
  challenge.alpha = hf_new(15);
  challenge.beta = hf_new(12);
  challenge.gamma = hf_new(13);
  challenge.z = hf_new(5);
  challenge.v = hf_new(12);

  // Benchmark proving
  start = get_time_ms();
  plonk_prove(&plonk, &constraints, &assignments, &challenge, rand);
  end = get_time_ms();
  printf("Large PLONK prove: %.2f ms\n", end - start);

  // Run multiple iterations to get better timing
  printf("\nRunning 5 iterations of large PLONK prove...\n");
  double total_time = 0.0;
  for (int i = 0; i < 5; i++) {
    start = get_time_ms();
    plonk_prove(&plonk, &constraints, &assignments, &challenge, rand);
    end = get_time_ms();
    total_time += (end - start);
    printf("  Iteration %d: %.2f ms\n", i+1, end - start);
  }
  printf("Average large PLONK prove time: %.2f ms\n", total_time / 5.0);

  // Clean up
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
  printf("=== Large PLONK Benchmark Results ===\n");
  benchmark_large_plonk();
  printf("==============================\n");
  return 0;
}