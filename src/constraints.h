#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hf.h"

// the gate enforces the constraint:
// q_l*a + q_r*b + q_o*c + q_m*a*b + q_c = 0
typedef struct {
  HF q_l;
  HF q_r;
  HF q_o;
  HF q_m;
  HF q_c;
} GATE;

typedef enum {
  COPYOF_A,
  COPYOF_B,
  COPYOF_C
} COPY_OF_TYPE;

// represents references to wire assignments in the circuit
// (inputs and outputs)
typedef struct {
  COPY_OF_TYPE type;
  size_t index;
} COPY_OF;

// holds vectors of gate coefficients and copy constraints.
// provides methods to evaluate expressions and check
// if a set of assignments satisfies the constraints
typedef struct {
  HF *q_l;
  HF *q_r;
  HF *q_o;
  HF *q_m;
  HF *q_c;
  size_t num_gates;

  COPY_OF *c_a;
  COPY_OF *c_b;
  COPY_OF *c_c;
  size_t num_constraints;
} CONSTRAINTS;

// assignment and assignments structs represent individual and
// collections of wire assignments in the circuit.
typedef struct {
  HF a;
  HF b;
  HF c;
} ASSIGNMENT;

typedef struct {
  HF *a;
  HF *b;
  HF *c;
  size_t len;
} ASSIGNMENTS;

typedef enum {
  EXPR_VAR,
  EXPR_CONST,
  EXPR_SUM,
  EXPR_SUB,
  EXPR_MUL
} EXPR_TYPE;

typedef struct expression {
  EXPR_TYPE type;
  union {
    const char *var_name; // for EXPR_VAR
    HF const_value;    // for EXPR_CONST
    struct {
      struct expression *left;
      struct expression *right;
    } binary; // for EXPR_SUM, EXPR_SUB, EXPR_MUL
  } data;
} EXPRESSION;

GATE gate_new(HF q_l, HF q_r, HF q_o, HF q_m, HF q_c) {
  GATE g = {q_l, q_r, q_o, q_m, q_c};
  return g;
}

// sum gate: a + b - c = 0
GATE gate_sum_a_b() {
  return gate_new(hf_one(), hf_one(), hf_neg(hf_one()), hf_zero(), hf_zero());
}

// sub gate: a - b - c = 0
GATE gate_sub_a_b() {
  return gate_new(hf_one(), hf_neg(hf_one()), hf_neg(hf_one()), hf_zero(), hf_zero());
}

// mul gate: a * b - c = 0
GATE gate_mul_a_b() {
  return gate_new(hf_zero(), hf_zero(), hf_neg(hf_one()), hf_one(), hf_zero());
}

// binds a variable to a constant: a + q_c = 0
GATE gate_bind_a(HF value) {
  return gate_new(hf_one(), hf_zero(), hf_zero(), hf_zero(), value);
}

// Gate to bind variable to zero: c = 0
GATE gate_bind_to_zero() {
  return gate_new(hf_zero(), hf_zero(), hf_one(), hf_zero(), hf_zero());
}

CONSTRAINTS constraints_new(GATE *gates, size_t num_gates, COPY_OF *c_a, COPY_OF *c_b, COPY_OF *c_c, size_t num_constraints) {
  CONSTRAINTS cons;
  cons.q_l = (HF *)malloc(num_gates * sizeof(HF));
  cons.q_r = (HF *)malloc(num_gates * sizeof(HF));
  cons.q_o = (HF *)malloc(num_gates * sizeof(HF));
  cons.q_m = (HF *)malloc(num_gates * sizeof(HF));
  cons.q_c = (HF *)malloc(num_gates * sizeof(HF));
  cons.num_gates = num_gates;

  for (size_t i = 0; i < num_gates; i++) {
    cons.q_l[i] = gates[i].q_l;
    cons.q_r[i] = gates[i].q_r;
    cons.q_o[i] = gates[i].q_o;
    cons.q_m[i] = gates[i].q_m;
    cons.q_c[i] = gates[i].q_c;
  }

  cons.c_a = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));
  cons.c_b = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));
  cons.c_c = (COPY_OF *)malloc(num_constraints * sizeof(COPY_OF));
  cons.num_constraints = num_constraints;

  for (size_t i = 0; i < num_constraints; i++) {
    cons.c_a[i] = c_a[i];
    cons.c_b[i] = c_b[i];
    cons.c_c[i] = c_c[i];
  }

  return cons;
}

bool constraints_satisfy(const CONSTRAINTS *c, const ASSIGNMENTS *a) {
  size_t n = c->num_constraints;
  for (size_t i = 0; i < n; i ++) {
    // compute lhs = (q_l[i] * a[i]) + (q_r[i] * b[i]) + (q_o[i] * c[i])
    //             + (q_m[i] * a[i] * b[i]) + q_c[i]
    HF lhs = hf_zero();
    HF term1 = hf_mul(c->q_l[i], a->a[i]);
    HF term2 = hf_mul(c->q_r[i], a->b[i]);
    HF term3 = hf_mul(c->q_o[i], a->c[i]);
    HF a_b = hf_mul(a->a[i], a->b[i]);
    HF term4 = hf_mul(c->q_m[i], a_b);

    lhs = hf_add(lhs, term1);
    lhs = hf_add(lhs, term2);
    lhs = hf_add(lhs, term3);
    lhs = hf_add(lhs, term4);
    lhs = hf_add(lhs, c->q_c[i]);

    // check if lhs is zero
    if (!hf_equal(lhs, hf_zero())) {
      printf("Constraint %zu not satisfied.\n", i);
      return false;
    }
  }

  return true;
}

void constraints_free(CONSTRAINTS *cons) {
  free(cons->q_l);
  free(cons->q_r);
  free(cons->q_o);
  free(cons->q_m);
  free(cons->q_c);

  free(cons->c_a);
  free(cons->c_b);
  free(cons->c_c);
}

#define MAX_VARS 100

typedef struct {
  char *names[MAX_VARS];
  size_t indices[MAX_VARS];
  size_t count;
} VAR_MAP;

void var_map_init(VAR_MAP *vm) {
  vm->count = 0;
}

size_t var_map_get_or_add(VAR_MAP *vm, const char *name) {
  for (size_t i = 0; i < vm->count; i++) {
    if (strcmp(vm->names[i], name) == 0)
      return vm->indices[i];
  }

  // Duplicate the name to allocate memory
  vm->names[vm->count] = strdup(name);
  if (!vm->names[vm->count]) {
    fprintf(stderr, "Memory allocation failed in var_map_get_or_add\n");
    exit(EXIT_FAILURE);
  }

  vm->indices[vm->count] = vm->count;
  vm->count++;
  return vm->indices[vm->count - 1];
}

const char *var_map_get_name(VAR_MAP *vm, size_t index) {
  if (index < vm->count)
    return vm->names[index];
  return NULL;
}

void var_map_free(VAR_MAP *vm) {
  for (size_t i = 0; i < vm->count; i++)
    free(vm->names[i]);
  vm->count = 0;
}

typedef struct {
  GATE *gates;
  size_t *a_indices;
  size_t *b_indices;
  size_t *c_indices;
  size_t num_gates;
  size_t capacity;
} GATE_LIST;

void gate_list_init(GATE_LIST *gl) {
  gl->capacity = 10;
  gl->gates = (GATE *)malloc(gl->capacity * sizeof(GATE));
  gl->a_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->b_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->c_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->num_gates = 0;
}

void gate_list_append(GATE_LIST *gl, GATE g, size_t a_index, size_t b_index, size_t c_index) {
  if (gl->num_gates >= gl->capacity) {
    gl->capacity *= 2;
    gl->gates = (GATE *)realloc(gl->gates, gl->capacity * sizeof(GATE));
    gl->a_indices = (size_t *)realloc(gl->a_indices, gl->capacity * sizeof(size_t));
    gl->b_indices = (size_t *)realloc(gl->b_indices, gl->capacity * sizeof(size_t));
    gl->c_indices = (size_t *)realloc(gl->c_indices, gl->capacity * sizeof(size_t));

    // Check for allocation failure
    if (!gl->gates || !gl->a_indices || !gl->b_indices || !gl->c_indices) {
      fprintf(stderr, "Memory reallocation failed in gate_list_append\n");
      exit(EXIT_FAILURE);
    }
  }
  gl->gates[gl->num_gates] = g;
  gl->a_indices[gl->num_gates] = a_index;
  gl->b_indices[gl->num_gates] = b_index;
  gl->c_indices[gl->num_gates] = c_index;
  gl->num_gates++;
}

void gate_list_free(GATE_LIST *gl) {
  free(gl->gates);
  free(gl->a_indices);
  free(gl->b_indices);
  free(gl->c_indices);
}

size_t eval_expr(EXPRESSION *expr, VAR_MAP *vars, GATE_LIST *gates) {
  switch (expr->type) {
  case EXPR_VAR: {
    return var_map_get_or_add(vars, expr->data.var_name);
  }
  case EXPR_CONST: {
    // constants need to be handled appropriately; for simplicity, we can assign them as new variables
    char const_name[20];
    snprintf(const_name, sizeof(const_name), "const_%u", expr->data.const_value.value);
    return var_map_get_or_add(vars, const_name);
  }
  case EXPR_SUM:
  case EXPR_SUB:
  case EXPR_MUL: {
    size_t l = eval_expr(expr->data.binary.left, vars, gates);
    size_t r = eval_expr(expr->data.binary.right, vars, gates);
    size_t n = vars->count;
    char var_name[20];
    snprintf(var_name, sizeof(var_name), "v%zu", n);
    var_map_get_or_add(vars, var_name);

    GATE gate;
    if (expr->type == EXPR_SUM) {
      gate = gate_sum_a_b();
    } else if (expr->type == EXPR_SUB) {
      gate = gate_sub_a_b();
    } else {
      gate = gate_mul_a_b();
    }
    gate_list_append(gates, gate, l, r, n);
    return n;
  }
  default:
    fprintf(stderr, "Unknown expression type\n");
    exit(EXIT_FAILURE);
  }
}

#endif // CONSTRAINTS_H
