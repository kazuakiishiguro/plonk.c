#ifndef PLONK_H
#define PLONK_H

#include <stdio.h>
#include <stdlib.h>

#include "pairing.h"
#include "poly.h"

typedef struct {
  // fields
  u8_fe gf;
  u8_fe hf;

  // ec points
  g1_p g1_generator;
  g2_p g2_generator;

  // constants
  u8_fe k1; // k1 x omega coset generator
  u8_fe k2; // k2 x omega coset generator
  u8_fe omega; // the generator in hf

  // function to map hf to gf
  u8_fe (*g_f)(u8_fe);
} plonk_types;

#endif // PLONK_H
