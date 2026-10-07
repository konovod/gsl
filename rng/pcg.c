/* rng/pcg.c
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

/* PCG32, the 32-bit output member of the PCG family:
 *
 *   M. E. O'Neill, "PCG: A Family of Simple Fast Space-Efficient
 *   Statistically Good Algorithms for Random Number Generation",
 *   Harvey Mudd College, 2014.
 *
 * The reference implementation is distributed under the Apache License
 * 2.0, which is compatible with the GPL; this is an independent C
 * implementation of the published algorithm.
 *
 * The state is a 64-bit LCG advanced by
 *
 *   state <- state * 6364136223846793005 + inc
 *
 * with an odd increment inc, followed by the output permutation
 *
 *   xorshifted = ((state >> 18) ^ state) >> 27
 *   rot        = state >> 59
 *   output     = (xorshifted >> rot) | (xorshifted << ((-rot) & 31))
 *
 * The period is 2^64.
 */

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>

static inline unsigned long int pcg32_get (void *vstate);
static double pcg32_get_double (void *vstate);
static void pcg32_set (void *vstate, unsigned long int s);

typedef struct
{
  uint64_t state;
  uint64_t inc;
}
pcg32_state_t;

static inline unsigned long int
pcg32_get (void *vstate)
{
  pcg32_state_t * state = (pcg32_state_t *) vstate;
  uint64_t oldstate = state->state;
  uint32_t xorshifted, rot;

  state->state = oldstate * UINT64_C(6364136223846793005) + state->inc;

  xorshifted = (uint32_t) (((oldstate >> 18) ^ oldstate) >> 27);
  rot = (uint32_t) (oldstate >> 59);

  return (unsigned long int)
    ((xorshifted >> rot) | (xorshifted << ((-rot) & 31)));
}

static double
pcg32_get_double (void *vstate)
{
  return pcg32_get (vstate) / 4294967296.0;
}

static void
pcg32_set (void *vstate, unsigned long int s)
{
  pcg32_state_t * state = (pcg32_state_t *) vstate;

  /* Seed both the LCG state and the stream selector from s, following
     the reference pcg32_srandom_r(state, s, s) so that distinct seeds
     give distinct streams. */

  state->state = 0;
  state->inc = ((uint64_t) s << 1) | 1;
  pcg32_get (state);
  state->state += (uint64_t) s;
  pcg32_get (state);
}

static const gsl_rng_type pcg32_type =
{"pcg32",                       /* name */
 0xffffffffUL,                  /* RAND_MAX */
 0,                             /* RAND_MIN */
 sizeof (pcg32_state_t),
 &pcg32_set,
 &pcg32_get,
 &pcg32_get_double};

const gsl_rng_type *gsl_rng_pcg32 = &pcg32_type;
