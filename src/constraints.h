#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fe.h"

// the gate enforces the constraint:
// q_l*a + q_r*b + q_o*c + q_m*a*b + q_c = 0
typedef struct {
  u8_fe q_l;
  u8_fe q_r;
  u8_fe q_o;
  u8_fe q_m;
  u8_fe q_c;
} gate;

typedef enum {
  COPYOF_A,
  COPYOF_B,
  COPYOF_C
} copy_of_type;

// represents references to wire assignments in the circuit
// (inputs and outputs)
typedef struct {
  copy_of_type type;
  size_t index;
} copy_of;

// holds vectors of gate coefficients and copy constraints.
// provides methods to evaluate expressions and check
// if a set of assignments satisfies the constraints
typedef struct {
  u8_fe *q_l;
  u8_fe *q_r;
  u8_fe *q_o;
  u8_fe *q_m;
  u8_fe *q_c;
  size_t num_gates;

  copy_of *c_a;
  copy_of *c_b;
  copy_of *c_c;
  size_t num_constraints;
} constraints;

// assignment and assignments structs represent individual and
// collections of wire assignments in the circuit.
typedef struct {
  u8_fe a;
  u8_fe b;
  u8_fe c;
} assignment;

typedef struct {
  u8_fe *a;
  u8_fe *b;
  u8_fe *c;
  size_t len;
} assignments;

typedef enum {
  EXPR_VAR,
  EXPR_CONST,
  EXPR_SUM,
  EXPR_SUB,
  EXPR_MUL
} expr_type;

typedef struct expression {
  expr_type type;
  union {
    const char *var_name; // for EXPR_VAR
    u8_fe const_value;    // for EXPR_CONST
    struct {
      struct expression *left;
      struct expression *right;
    } binary; // for EXPR_SUM, EXPR_SUB, EXPR_MUL
  } data;
} expression;

gate gate_new(u8_fe q_l, u8_fe q_r, u8_fe q_o, u8_fe q_m, u8_fe q_c) {
  gate g = {q_l, q_r, q_o, q_m, q_c};
  return g;
}

// sum gate: a + b - c = 0
gate gate_sum_a_b() {
  return gate_new(u8_fe_one(), u8_fe_one(), u8_fe_neg(u8_fe_one()), u8_fe_zero(), u8_fe_zero());
}

// sub gate: a - b - c = 0
gate gate_sub_a_b() {
  return gate_new(u8_fe_one(), u8_fe_neg(u8_fe_one()), u8_fe_neg(u8_fe_one()), u8_fe_zero(), u8_fe_zero());
}

// mul gate: a * b - c = 0
gate gate_mul_a_b() {
  return gate_new(u8_fe_zero(), u8_fe_zero(), u8_fe_neg(u8_fe_one()), u8_fe_one(), u8_fe_zero());
}

// binds a variable to a constant: a + q_c = 0
gate gate_bind_a(u8_fe value) {
  return gate_new(u8_fe_one(), u8_fe_zero(), u8_fe_zero(), u8_fe_zero(), value);
}

// Gate to bind variable to zero: c = 0
gate gate_bind_to_zero() {
  return gate_new(u8_fe_zero(), u8_fe_zero(), u8_fe_one(), u8_fe_zero(), u8_fe_zero());
}

constraints constraints_new(gate *gates, size_t num_gates, copy_of *c_a, copy_of *c_b, copy_of *c_c, size_t num_constraints) {
  constraints cons;
  cons.q_l = (u8_fe *)malloc(num_gates * sizeof(u8_fe));
  cons.q_r = (u8_fe *)malloc(num_gates * sizeof(u8_fe));
  cons.q_o = (u8_fe *)malloc(num_gates * sizeof(u8_fe));
  cons.q_m = (u8_fe *)malloc(num_gates * sizeof(u8_fe));
  cons.q_c = (u8_fe *)malloc(num_gates * sizeof(u8_fe));
  cons.num_gates = num_gates;

  for (size_t i = 0; i < num_gates; i++) {
    cons.q_l[i] = gates[i].q_l;
    cons.q_r[i] = gates[i].q_r;
    cons.q_o[i] = gates[i].q_o;
    cons.q_m[i] = gates[i].q_m;
    cons.q_c[i] = gates[i].q_c;
  }

  cons.c_a = (copy_of *)malloc(num_constraints * sizeof(copy_of));
  cons.c_b = (copy_of *)malloc(num_constraints * sizeof(copy_of));
  cons.c_c = (copy_of *)malloc(num_constraints * sizeof(copy_of));
  cons.num_constraints = num_constraints;

  for (size_t i = 0; i < num_constraints; i++) {
    cons.c_a[i] = c_a[i];
    cons.c_b[i] = c_b[i];
    cons.c_c[i] = c_c[i];
  }

  return cons;
}


bool constraints_satisfy(const constraints *c, const assignments *a) {
  size_t n = c->num_constraints;
  for (size_t i = 0; i < n; i ++) {
    // compute lhs = (q_l[i] * a[i]) + (q_r[i] * b[i]) + (q_o[i] * c[i])
    //             + (q_m[i] * a[i] * b[i]) + q_c[i]
    u8_fe lhs = u8_fe_zero();
    u8_fe term1 = u8_fe_mul(c->q_l[i], a->a[i]);
    u8_fe term2 = u8_fe_mul(c->q_r[i], a->b[i]);
    u8_fe term3 = u8_fe_mul(c->q_o[i], a->c[i]);
    u8_fe a_b = u8_fe_mul(a->a[i], a->b[i]);
    u8_fe term4 = u8_fe_mul(c->q_m[i], a_b);

    lhs = u8_fe_add(lhs, term1);
    lhs = u8_fe_add(lhs, term2);
    lhs = u8_fe_add(lhs, term3);
    lhs = u8_fe_add(lhs, term4);
    lhs = u8_fe_add(lhs, c->q_c[i]);

    // check if lhs is zero
    if (!u8_fe_equal(lhs, u8_fe_zero())) {
      return false;
    }
  }

  return true;
}

void constraints_free(constraints *cons) {
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
} var_map;

void var_map_init(var_map *vm) {
  vm->count = 0;
}

size_t var_map_get_or_add(var_map *vm, const char *name) {
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

const char *var_map_get_name(var_map *vm, size_t index) {
  if (index < vm->count)
    return vm->names[index];
  return NULL;
}

void var_map_free(var_map *vm) {
  for (size_t i = 0; i < vm->count; i++)
    free(vm->names[i]);
  vm->count = 0;
}

typedef struct {
  gate *gates;
  size_t *a_indices;
  size_t *b_indices;
  size_t *c_indices;
  size_t num_gates;
  size_t capacity;
} gate_list;

void gate_list_init(gate_list *gl) {
  gl->capacity = 10;
  gl->gates = (gate *)malloc(gl->capacity * sizeof(gate));
  gl->a_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->b_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->c_indices = (size_t *)malloc(gl->capacity * sizeof(size_t));
  gl->num_gates = 0;
}

void gate_list_append(gate_list *gl, gate g, size_t a_index, size_t b_index, size_t c_index) {
  if (gl->num_gates >= gl->capacity) {
    gl->capacity *= 2;
    gl->gates = (gate *)realloc(gl->gates, gl->capacity * sizeof(gate));
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

void gate_list_free(gate_list *gl) {
  free(gl->gates);
  free(gl->a_indices);
  free(gl->b_indices);
  free(gl->c_indices);
}

size_t eval_expr(expression *expr, var_map *vars, gate_list *gates) {
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

    gate g;
    if (expr->type == EXPR_SUM) {
      g = gate_sum_a_b();
    } else if (expr->type == EXPR_SUB) {
      g = gate_sub_a_b();
    } else {
      g = gate_mul_a_b();
    }
    gate_list_append(gates, g, l, r, n);
    return n;
  }
  default:
    fprintf(stderr, "Unknown expression type\n");
    exit(EXIT_FAILURE);
  }
}

#endif // CONSTRAINTS_H
