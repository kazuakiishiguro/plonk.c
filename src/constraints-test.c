#include <assert.h>

#include "constraints.h"

#define UNUSED_INDEX ((size_t)(-1))

void test_expr() {
  // initialize variable map and gate list
  var_map vars;
  var_map_init(&vars);

  gate_list gates;
  gate_list_init(&gates);

  // define expressions
  expression *a = (expression *)malloc(sizeof(expression));
  a->type = EXPR_VAR;
  a->data.var_name = "a";

  expression *b = (expression *)malloc(sizeof(expression));
  b->type = EXPR_VAR;
  b->data.var_name = "b";

  expression *c = (expression *)malloc(sizeof(expression));
  c->type = EXPR_VAR;
  c->data.var_name = "c";

  // construct the expression: (a * a) + (b * b) - (c * c)
  expression *aa = (expression *)malloc(sizeof(expression));
  aa->type = EXPR_MUL;
  aa->data.binary.left = a;
  aa->data.binary.right = a;

  expression *bb = (expression *)malloc(sizeof(expression));
  bb->type = EXPR_MUL;
  bb->data.binary.left = b;
  bb->data.binary.right = b;

  expression *cc = (expression *)malloc(sizeof(expression));
  cc->type = EXPR_MUL;
  cc->data.binary.left = c;
  cc->data.binary.right = c;

  expression *sum = (expression *)malloc(sizeof(expression));
  sum->type = EXPR_SUM;
  sum->data.binary.left = aa;
  sum->data.binary.right = bb;

  expression *pitagoras = (expression *)malloc(sizeof(expression));
  pitagoras->type = EXPR_SUB;
  pitagoras->data.binary.left = sum;
  pitagoras->data.binary.right = cc;

  // evaluate the expression and generate gates
  size_t result_index = eval_expr(pitagoras, &vars, &gates);

  // add a gate to bind the result to zero
  gate bind_gate = gate_bind_to_zero();
  gate_list_append(&gates, bind_gate, UNUSED_INDEX, UNUSED_INDEX, result_index);

  // print variable mappings
  printf("Variable mappings:\n");
  for (size_t i = 0; i < vars.count; i++)
    printf("%zu => %s\n", i, vars.names[i]);

  // print gates
  printf("\nGenerated gates:\n");
  for (size_t i = 0; i < gates.num_gates; i++) {
    gate g = gates.gates[i];
    size_t a_idx = gates.a_indices[i];
    size_t b_idx = gates.b_indices[i];
    size_t c_idx = gates.c_indices[i];

    int q_l = (g.q_l.value <= MODULO / 2) ? g.q_l.value : g.q_l.value - MODULO;
    int q_r = (g.q_r.value <= MODULO / 2) ? g.q_r.value : g.q_r.value - MODULO;
    int q_o = (g.q_o.value <= MODULO / 2) ? g.q_o.value : g.q_o.value - MODULO;
    int q_m = (g.q_m.value <= MODULO / 2) ? g.q_m.value : g.q_m.value - MODULO;
    int q_c = (g.q_c.value <= MODULO / 2) ? g.q_c.value : g.q_c.value - MODULO;

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

  constraints c;
  c.num_constraints = n;
  c.q_l = (u8_fe *)malloc(n * sizeof(u8_fe));
  c.q_r = (u8_fe *)malloc(n * sizeof(u8_fe));
  c.q_o = (u8_fe *)malloc(n * sizeof(u8_fe));
  c.q_m = (u8_fe *)malloc(n * sizeof(u8_fe));
  c.q_c = (u8_fe *)malloc(n * sizeof(u8_fe));

  if (!c.q_l || !c.q_r || !c.q_o || !c.q_m || !c.q_c) {
    fprintf(stderr, "Memory allocation failed in test_constraints_satisfy\n");
    exit(EXIT_FAILURE);
  }

  // set up the constraint coefficients for c = a + b
  // q_l = 1, q_r = 1, q_o = -1, q_m = 0, q_c = 0
  c.q_l[0] = u8_fe_one();
  c.q_r[0] = u8_fe_one();
  c.q_o[0] = u8_fe_neg(u8_fe_one()); // -1 mod MODULO
  c.q_m[0] = u8_fe_zero();
  c.q_c[0] = u8_fe_zero();

  // Assignments
  assignments a;
  a.len = n;
  a.a = (u8_fe *)malloc(n * sizeof(u8_fe));
  a.b = (u8_fe *)malloc(n * sizeof(u8_fe));
  a.c = (u8_fe *)malloc(n * sizeof(u8_fe));

  if (!a.a || !a.b || !a.c) {
    fprintf(stderr, "Memory allocation failed in test_constraints_satisfy\n");
    exit(EXIT_FAILURE);
  }

  // test case 1: correct assignments (a = 2, b = 3, c = 5)
  a.a[0] = u8_fe_new(2);
  a.b[0] = u8_fe_new(3);
  a.c[0] = u8_fe_new(5);

  // check if constraints are satisfied
  assert(constraints_satisfy(&c, &a));

  // test case 2: incorrect assignments (a = 2, b = 3, c = 6)
  a.c[0] = u8_fe_new(6);
  assert(!constraints_satisfy(&c, &a));

  constraints_free(&c);
  free(a.a);
  free(a.b);
  free(a.c);
}

int main() {
  test_expr();
  test_constraints_satisfy();
  return 0;
}
