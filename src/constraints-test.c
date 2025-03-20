#include <assert.h>

#include "constraints.h"

#define UNUSED_INDEX ((size_t)(-1))

void test_expr() {
  // initialize variable map and gate list
  VAR_MAP vars;
  var_map_init(&vars);

  GATE_LIST gates;
  gate_list_init(&gates);

  // define expressions
  EXPRESSION *a = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  a->type = EXPR_VAR;
  a->data.var_name = "a";

  EXPRESSION *b = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  b->type = EXPR_VAR;
  b->data.var_name = "b";

  EXPRESSION *c = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  c->type = EXPR_VAR;
  c->data.var_name = "c";

  // construct the expression: (a * a) + (b * b) - (c * c)
  EXPRESSION *aa = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  aa->type = EXPR_MUL;
  aa->data.binary.left = a;
  aa->data.binary.right = a;

  EXPRESSION *bb = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  bb->type = EXPR_MUL;
  bb->data.binary.left = b;
  bb->data.binary.right = b;

  EXPRESSION *cc = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  cc->type = EXPR_MUL;
  cc->data.binary.left = c;
  cc->data.binary.right = c;

  EXPRESSION *sum = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  sum->type = EXPR_SUM;
  sum->data.binary.left = aa;
  sum->data.binary.right = bb;

  EXPRESSION *pitagoras = (EXPRESSION *)malloc(sizeof(EXPRESSION));
  pitagoras->type = EXPR_SUB;
  pitagoras->data.binary.left = sum;
  pitagoras->data.binary.right = cc;

  // evaluate the EXPRESSION and generate gates
  size_t result_index = eval_expr(pitagoras, &vars, &gates);

  // add a gate to bind the result to zero
  GATE bind_gate = gate_bind_to_zero();
  gate_list_append(&gates, bind_gate, UNUSED_INDEX, UNUSED_INDEX, result_index);

  // print variable mappings
  printf("Variable mappings:\n");
  for (size_t i = 0; i < vars.count; i++)
    printf("%zu => %s\n", i, vars.names[i]);

  // print gates
  printf("\nGenerated gates:\n");
  for (size_t i = 0; i < gates.num_gates; i++) {
    GATE gate = gates.gates[i];
    size_t a_idx = gates.a_indices[i];
    size_t b_idx = gates.b_indices[i];
    size_t c_idx = gates.c_indices[i];

    int q_l = (gate.q_l.value <= MODULO_HF / 2) ? gate.q_l.value : gate.q_l.value - MODULO_HF;
    int q_r = (gate.q_r.value <= MODULO_HF / 2) ? gate.q_r.value : gate.q_r.value - MODULO_HF;
    int q_o = (gate.q_o.value <= MODULO_HF / 2) ? gate.q_o.value : gate.q_o.value - MODULO_HF;
    int q_m = (gate.q_m.value <= MODULO_HF / 2) ? gate.q_m.value : gate.q_m.value - MODULO_HF;
    int q_c = (gate.q_c.value <= MODULO_HF / 2) ? gate.q_c.value : gate.q_c.value - MODULO_HF;

    const char *a_name = (a_idx != UNUSED_INDEX) ? vars.names[a_idx] : "(unused)";
    const char *b_name = (b_idx != UNUSED_INDEX) ? vars.names[b_idx] : "(unused)";
    const char *c_name = vars.names[c_idx];

    printf("Gate %zu: %d*a + %d*b + %d*c + %d*a*b + %d = 0\n",
	   i, q_l, q_r, q_o, q_m, q_c);

    printf("Wires: a = %s, b = %s, c = %s\n",
	   a_name,
	   b_name,
	   c_name);
  }

  printf("\n");

  gate_list_free(&gates);
  var_map_free(&vars);

  // Free expressions (avoid memory leaks)
  free(pitagoras); // pitagoras depends on sum and cc
  free(sum);       // sum depends on aa and bb
  free(aa);        // aa depends on a and a
  free(bb);        // bb depends on b and b
  free(cc);        // cc depends on c and c
  free(a);         // a
  free(b);         // b
  free(c);         // c
}

void test_constraints_satisfy() {
  // simple constraint system representing c = a + b
  size_t n = 1; // num of constraints

  CONSTRAINTS constraints;
  constraints.num_constraints = n;
  constraints.q_l = (HF *)malloc(n * sizeof(HF));
  constraints.q_r = (HF *)malloc(n * sizeof(HF));
  constraints.q_o = (HF *)malloc(n * sizeof(HF));
  constraints.q_m = (HF *)malloc(n * sizeof(HF));
  constraints.q_c = (HF *)malloc(n * sizeof(HF));

  if (!constraints.q_l || !constraints.q_r || !constraints.q_o || !constraints.q_m || !constraints.q_c) {
    fprintf(stderr, "Memory allocation failed in test_constraints_satisfy\n");
    exit(EXIT_FAILURE);
  }

  // set up the constraint coefficients for c = a + b
  // q_l = 1, q_r = 1, q_o = -1, q_m = 0, q_c = 0
  constraints.q_l[0] = hf_one();
  constraints.q_r[0] = hf_one();
  constraints.q_o[0] = hf_neg(hf_one()); // -1 mod MODULO_HF
  constraints.q_m[0] = hf_zero();
  constraints.q_c[0] = hf_zero();

  // Assignments
  ASSIGNMENTS assignments;
  assignments.len = n;
  assignments.a = (HF *)malloc(n * sizeof(HF));
  assignments.b = (HF *)malloc(n * sizeof(HF));
  assignments.c = (HF *)malloc(n * sizeof(HF));

  if (!assignments.a || !assignments.b || !assignments.c) {
    fprintf(stderr, "Memory allocation failed in test_constraints_satisfy\n");
    exit(EXIT_FAILURE);
  }

  // test case 1: correct assignments (a = 2, b = 3, c = 5)
  assignments.a[0] = hf_new(2);
  assignments.b[0] = hf_new(3);
  assignments.c[0] = hf_new(5);

  // check if constraints are satisfied
  assert(constraints_satisfy(&constraints, &assignments));

  // test case 2: incorrect assignments (a = 2, b = 3, c = 6)
  assignments.c[0] = hf_new(6);
  assert(!constraints_satisfy(&constraints, &assignments));

  free(assignments.a);
  free(assignments.b);
  free(assignments.c);
}

int main() {
  test_expr();
  test_constraints_satisfy();
  return 0;
}
