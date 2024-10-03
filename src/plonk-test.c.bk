#include <assert.h>
#include "plonk.h"

/* void test_plonk() { */
/*   // step 1: create trustet setup (srs) */
/*   // toxic waste */
/*   u8_fe secret = u8_fe_new(2); */
/*   // number of points in srs */
/*   size_t n = 10; */

/*   srs s = srs_create(secret, n); */

/*   // step 2: create the plonk instance */
/*   // number of omega powers (size of H) */
/*   size_t num_constraints = 4; */
/*   size_t omega_pows = num_constraints + 1; */
/*   plonk p = plonk_new(s, omega_pows); */

/*   // step 3: set up constraints and assignments */
/*   // define the gates */
/*   // gates: */
/*   // - gate 1: mul a * b */
/*   // - gate 2: mul a * b */
/*   // - gate 3: mul a * b */
/*   // - gate 4: mul a + b */

/*   constraints c; */
/*   c.num_constraints = num_constraints; */

/*   // initialize q_m, q_l, q_r, q_o, q_c arrays */
/*   c.q_m = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   c.q_l = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   c.q_r = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   c.q_o = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   c.q_c = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */

/*   if (!c.q_m || !c.q_l || !c.q_r || !c.q_o || !c.q_c) { */
/*     fprintf(stderr, "Memory allocation failed in test_plonk\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */

/*   // date definitions (for simplicity, we can define the gates manually) */
/*   // for multiplication gates: q_m = l, q_o = -l */
/*   // for addition gate: q_l = l, q_r = l, q_o = -l */

/*   // gate 1: mul a * b */
/*   c.q_m[0] = u8_fe_new(1); */
/*   c.q_l[0] = u8_fe_new(0); */
/*   c.q_r[0] = u8_fe_new(0); */
/*   c.q_o[0] = u8_fe_neg(u8_fe_new(1)); */
/*   c.q_c[0] = u8_fe_new(0); */

/*   // gate 2: mul a * b */
/*   c.q_m[1] = u8_fe_new(1); */
/*   c.q_l[1] = u8_fe_new(0); */
/*   c.q_r[1] = u8_fe_new(0); */
/*   c.q_o[1] = u8_fe_neg(u8_fe_new(1)); */
/*   c.q_c[1] = u8_fe_new(0); */

/*   // gate 3: mul a * b */
/*   c.q_m[2] = u8_fe_new(1); */
/*   c.q_l[2] = u8_fe_new(0); */
/*   c.q_r[2] = u8_fe_new(0); */
/*   c.q_o[2] = u8_fe_neg(u8_fe_new(1)); */
/*   c.q_c[2] = u8_fe_new(0); */

/*   // gate 4: sum a + b */
/*   c.q_m[3] = u8_fe_new(0); */
/*   c.q_l[3] = u8_fe_new(1); */
/*   c.q_r[3] = u8_fe_new(1); */
/*   c.q_o[3] = u8_fe_neg(u8_fe_new(1)); */
/*   c.q_c[3] = u8_fe_new(0); */

/*   // copy constraints (permutation) */
/*   // for simplicity, let's assume c_a, c_b, c_c are arrays of copy_of structures */
/*   c.c_a = (copy_of *)malloc(num_constraints * sizeof(copy_of)); */
/*   c.c_b = (copy_of *)malloc(num_constraints * sizeof(copy_of)); */
/*   c.c_c = (copy_of *)malloc(num_constraints * sizeof(copy_of)); */

/*   if (!c.c_a || !c.c_b || !c.c_c) { */
/*     fprintf(stderr, "Memory allocation failed in test_plonk\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */

/*   // define copy constraints */
/*   // c_a: [b(1), b(2), b(3), c(1)] */
/*   c.c_a[0] = (copy_of){COPYOF_B, 1}; */
/*   c.c_a[1] = (copy_of){COPYOF_B, 2}; */
/*   c.c_a[2] = (copy_of){COPYOF_B, 3}; */
/*   c.c_a[3] = (copy_of){COPYOF_C, 1}; */

/*   // c_b: [a(1), a(2), a(3), c(2)] */
/*   c.c_b[0] = (copy_of){COPYOF_A, 1}; */
/*   c.c_b[1] = (copy_of){COPYOF_A, 2}; */
/*   c.c_b[2] = (copy_of){COPYOF_A, 3}; */
/*   c.c_b[3] = (copy_of){COPYOF_C, 2}; */

/*   // c_c: [a(4), a(4), a(4), c(3)] */
/*   c.c_c[0] = (copy_of){COPYOF_A, 4}; */
/*   c.c_c[1] = (copy_of){COPYOF_A, 4}; */
/*   c.c_c[2] = (copy_of){COPYOF_A, 4}; */
/*   c.c_c[3] = (copy_of){COPYOF_C, 3}; */

/*   // step 4: define assignments */
/*   assignments a; */
/*   a.len = num_constraints; */
/*   a.a = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   a.b = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */
/*   a.c = (u8_fe *)malloc(num_constraints * sizeof(u8_fe)); */

/*   if (!a.a || !a.b || !a.c) { */
/*     fprintf(stderr, "Memory allocation failed in test_plonk\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */

/*   // assign values */
/*   // assignment 1: a = 3, b = 3, c = 9 */
/*   a.a[0] = u8_fe_new(3); */
/*   a.b[0] = u8_fe_new(3); */
/*   a.c[0] = u8_fe_new(9); */

/*   // assignment 2: a = 4, b = 4, c = 16 */
/*   a.a[1] = u8_fe_new(4); */
/*   a.b[1] = u8_fe_new(4); */
/*   a.c[1] = u8_fe_new(16); */

/*   // assignment 3: a = 5, b = 5, c = 25 */
/*   a.a[2] = u8_fe_new(5); */
/*   a.b[2] = u8_fe_new(5); */
/*   a.c[2] = u8_fe_new(25); */

/*   // assignment 4: a = 9, b = 16, c = 25 */
/*   a.a[3] = u8_fe_new(9); */
/*   a.b[3] = u8_fe_new(16); */
/*   a.c[3] = u8_fe_new(25); */

/*   // step 5: define random numbers (the b's) */
/*   u8_fe rand[9] = { */
/*     u8_fe_new(7), */
/*     u8_fe_new(4), */
/*     u8_fe_new(11), */
/*     u8_fe_new(12), */
/*     u8_fe_new(16), */
/*     u8_fe_new(2), */
/*     u8_fe_new(14), */
/*     u8_fe_new(11), */
/*     u8_fe_new(7) */
/*   }; */

/*   // step 6: define challenge calues */
/*   challenge ch; */
/*   ch.alpha = u8_fe_new(15); */
/*   ch.beta = u8_fe_new(12); */
/*   ch.gamma = u8_fe_new(13); */
/*   ch.z = u8_fe_new(5); */
/*   ch.z = u8_fe_new(12); */

/*   // step 7: call plonk_prove to generate a proof */
/*   proof prf = plonk_prove(&p, &c, &a, &ch, rand); */

/*   // step 8: define the expected proof (for comparison) */
/*   proof expected; */
/*   expected.a_s = g1_p_new(u8_fe_new(91).value, u8_fe_new(66).value); */
/*   expected.b_s = g1_p_new(u8_fe_new(26).value, u8_fe_new(45).value); */
/*   expected.c_s = g1_p_new(u8_fe_new(91).value, u8_fe_new(35).value); */
/*   expected.z_s = g1_p_new(u8_fe_new(32).value, u8_fe_new(59).value); */
/*   expected.t_lo_s = g1_p_new(u8_fe_new(12).value, u8_fe_new(32).value); */
/*   expected.t_mid_s = g1_p_new(u8_fe_new(26).value, u8_fe_new(45).value); */
/*   expected.t_hi_s = g1_p_new(u8_fe_new(91).value, u8_fe_new(66).value); */
/*   expected.w_z_s = g1_p_new(u8_fe_new(91).value, u8_fe_new(35).value); */
/*   expected.w_z_omega_s = g1_p_new(u8_fe_new(65).value, u8_fe_new(98).value); */
/*   expected.a_z = u8_fe_new(15); */
/*   expected.b_z = u8_fe_new(13); */
/*   expected.c_z = u8_fe_new(5); */
/*   expected.s_sigma_1_z = u8_fe_new(1); */
/*   expected.s_sigma_2_z = u8_fe_new(12); */
/*   expected.r_z = u8_fe_new(15); */
/*   expected.z_omega_z = u8_fe_new(15); */

/*   // step 9: compare the generated proof with the expected proof */
/*   // TODO: check proof equality */

/*   // step 10: verify the proof */
/*   u8_fe rand_verifier[1] = { u8_fe_new(4) }; */
/*   assert(plonk_verify(&p, &c, &prf, &ch, rand_verifier)); */

/*   free(c.q_m); */
/*   free(c.q_l); */
/*   free(c.q_r); */
/*   free(c.q_o); */
/*   free(c.q_c); */
/*   free(c.c_a); */
/*   free(c.c_b); */
/*   free(c.c_c); */
/*   free(a.a); */
/*   free(a.b); */
/*   free(a.c); */
/*   plonk_free(&p); */
/* } */

int main() {
  /* test_plonk(); */
  return 0;
}
