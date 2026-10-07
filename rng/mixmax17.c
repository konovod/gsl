/* rng/mixmax17.c
 * 
 * Copyright (C) 2026 The GSL Team
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* MIXMAX with N = 17 over the Mersenne prime field GF(2^61 - 1).
 *
 * This is the variant chosen as the default by CLHEP and hence Geant4:
 *
 *   K. Savvidy, "The MIXMAX random number generator",
 *   Computer Physics Communications 196 (2015) 161-165;
 *
 *   K. Savvidy and G. Savvidy, "Spectrum and entropy of C-systems.
 *   MIXMAX random number generator",
 *   Chaos, Solitons & Fractals 91 (2016) 33-38;
 *
 *   G. K. Savvidy and N. G. Ter-Arutyunyan-Savvidy, "On the Monte Carlo
 *   simulation of physical systems", J. Comput. Phys. 97 (1991) 566-572.
 *
 * The state is a vector V[0..16] of residues modulo M61 = 2^61 - 1
 * together with their sum.  One iteration multiplies the vector by the
 * one-parameter MIXMAX matrix,
 *
 *   V[i] <- V[i] + sum_{j<i} V[j] * m^(i-j)   (mod M61),
 *
 * whose special structure lets a whole iteration run in O(N) with the
 * MULWU partial-sum trick below.  V[0] holds the checksum and is never
 * returned, so each iteration emits the N-1 values V[1..16]; the period
 * of the matrix is about 10^294.  The reference C implementation (LGPLv3)
 * is distributed with CLHEP and ROOT; this is an independent
 * implementation of the published algorithm.
 *
 * The native output is a 61-bit integer.  GSL's integer interface returns
 * `unsigned long int`, which is 32 bits on LLP64 as well as LP64, so
 * gsl_rng_get returns the low 32 bits of each draw -- the same reduction
 * CLHEP applies in its `operator unsigned int`.  gsl_rng_uniform uses the
 * full 61-bit value.  Seeding follows the reference `seed_spbox`.
 */

#include <config.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>

#define MIXMAX17_N          17
#define MIXMAX17_BITS       61
#define MIXMAX17_M61        UINT64_C(2305843009213693951)  /* 2^61 - 1 */
#define MIXMAX17_SPECIALMUL 36
#define MIXMAX17_MULT64     UINT64_C(6364136223846793005)

static inline unsigned long int mixmax17_get (void *vstate);
static double mixmax17_get_double (void *vstate);
static void mixmax17_set (void *vstate, unsigned long int s);

typedef struct
{
  uint64_t v[MIXMAX17_N];
  uint64_t sumtot;
  int counter;
}
mixmax17_state_t;

/* Reduce k modulo 2^61 - 1 in the Mersenne style.  For the sums the
   iterate loop forms (k < 3 * 2^61) the result is congruent to k and is
   at most 2^61; mixmax17_next trims the values it returns to the range
   [0, 2^61 - 1]. */
static inline uint64_t
mixmax17_mod (uint64_t k)
{
  return (k & MIXMAX17_M61) + (k >> MIXMAX17_BITS);
}

/* The one-parameter matrix folded into the running partial sum. */
static inline uint64_t
mixmax17_mulwu (uint64_t k)
{
  return ((k << MIXMAX17_SPECIALMUL) & MIXMAX17_M61)
    ^ (k >> (MIXMAX17_BITS - MIXMAX17_SPECIALMUL));
}

static void
mixmax17_iterate (mixmax17_state_t * state)
{
  uint64_t sumtot, ovflow, tempP, tempV;
  int i;

  tempV = state->sumtot;
  state->v[0] = tempV;
  sumtot = tempV;
  ovflow = 0;
  tempP = 0;

  for (i = 1; i < MIXMAX17_N; i++)
    {
      uint64_t tempPO = mixmax17_mulwu (tempP);
      tempP = mixmax17_mod (tempP + state->v[i]);
      tempV = mixmax17_mod (tempV + tempP + tempPO);
      state->v[i] = tempV;
      sumtot += tempV;
      if (sumtot < tempV)       /* wrapped past 2^64 */
        ovflow++;
    }

  /* Every overflow adds 2^64 = 2^3 * 2^61 = 8 (mod 2^61 - 1). */
  state->sumtot = mixmax17_mod (mixmax17_mod (sumtot) + (ovflow << 3));
}

static uint64_t
mixmax17_next (mixmax17_state_t * state)
{
  if (state->counter < MIXMAX17_N)
    {
      return mixmax17_mod (state->v[state->counter++]);
    }
  else
    {
      mixmax17_iterate (state);
      state->counter = 2;
      return mixmax17_mod (state->v[1]);
    }
}

static unsigned long int
mixmax17_get (void *vstate)
{
  mixmax17_state_t * state = (mixmax17_state_t *) vstate;
  return (unsigned long int) (mixmax17_next (state) & UINT64_C(0xffffffff));
}

static double
mixmax17_get_double (void *vstate)
{
  mixmax17_state_t * state = (mixmax17_state_t *) vstate;
  return (double) mixmax17_next (state) / 2305843009213693952.0;  /* 2^61 */
}

static void
mixmax17_set (void *vstate, unsigned long int s)
{
  mixmax17_state_t * state = (mixmax17_state_t *) vstate;
  uint64_t l;
  int i;

  if (s == 0)
    s = 1;                      /* default seed is 1 */

  /* seed_spbox: a 64-bit LCG with a half-word swap. */
  l = (uint64_t) s;
  for (i = 0; i < MIXMAX17_N; i++)
    {
      l *= MIXMAX17_MULT64;
      l = (l << 32) ^ (l >> 32);
      state->v[i] = l & MIXMAX17_M61;
    }

  state->sumtot = 0;
  for (i = 0; i < MIXMAX17_N; i++)
    state->sumtot = mixmax17_mod (state->sumtot + state->v[i]);

  state->counter = MIXMAX17_N;  /* force an iteration on the first draw */
}

static const gsl_rng_type mixmax17_type =
{"mixmax17",                    /* name */
 0xffffffffUL,                  /* RAND_MAX */
 0,                             /* RAND_MIN */
 sizeof (mixmax17_state_t),
 &mixmax17_set,
 &mixmax17_get,
 &mixmax17_get_double};

const gsl_rng_type *gsl_rng_mixmax17 = &mixmax17_type;
