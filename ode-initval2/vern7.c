/* ode-initval2/vern7.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2024 Christian Krueger
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

/* Runge-Kutta 7(6), Verner
 *
 * J.H. Verner, Numerically optimal Runge--Kutta pairs 
 * with interpolants. Numerical Algorithms, 53, (2010) pp. 383--396
 *
 * https://www.sfu.ca/~jverner/
 * "even more efficient" from July 2024.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "odeiv_util.h"
#include "step_utils.c"

/* Vern7 constants */

static const double Abar[] = {
  .5163520172057869163393251056217968836723e-1,
  0,
  0,
  .2767172535461648728769641534539952501983,
  .3374175285287150670818592701488271741753,
  .1884488267810967803491085059046161195540,
  24.54134121634868026791753618430192161716,
  -68.81190284469011946382716084194838780382,
  44.41634281776488378396776021757684795437,
  0
};

static const double A[] = {
  .5089676583692947576073561095512200263213e-1,
  0,
  0,
  .2793777374763233901369432426263934138476,
  .3281330142746535239936396881369403928344,
  .2241721218186151033581794837350130009230,
  .7874574778015076584344903106189416715189,
  0,
  0,
  -.6700371172080291516839883360724104817561
};

static const double ah[] = {
  .69e-1,
  .118,
  .177,
  .501,
  .7737799115305331003715765296862487670813,
  .994,
  .998
};

static const double b21 = .69e-1;

static const double b3[] = {
  .1710144927536231884057971014492753623188e-1,
  .1008985507246376811594202898550724637681
};

static const double b4[] = {
  .4425e-1,
  0,
  .13275
};

static const double b5[] = {
  .7353445130709566216604424016087331226659,
  0,
  -2.830160657856937661591496696351623096811,
  2.595816144785981039931054294742889974145
};

static const double b6[] = {
  -12.21580485360407974005910916471598682362,
  0,
  48.82665485823736062335980699373053427134,
  -38.55615592319928364666616600329792491404,
  2.719085830096535863737044703969626233400
};

static const double b7[] = {
  108.8614188704176574066699618897203578466,
  0,
  -432.4521181775777896358931629332707752654,
  343.9115281800118289547200158889409233641,
  -20.55041135925273709189369488701721016265,
  1.223582486401040366396880041626704217305
};

static const double b8[] = {
  113.4755131883738522204615568160304033854,
  0,
  -450.8122021555997002820400438087344405365,
  358.5132765190089889943579090008312808216,
  -21.45046667648445540174055882443151176550,
  1.274053318605952891766776667539031508649,
  -.2174193904638422805639851234763413667602e-2
};

static const double b9[] = {
  115.6996223324232534824963925993127275021,
  0,
  -459.6635446100248030478961869239726305957,
  365.5534717131745930309149378867953890507,
  -21.88511586349784824146225495848432937529,
  1.298718109698721459187976480852777474315,
  -.5318700918481883515898878747322241917739e-4,
  -.3098494764731864405706095716460833640254e-2
};

static const double b10[] = {
  124.1543935612464600014576130437603883332,
  0,
  -493.2318713314597046194663569971348299332,
  392.2086219315800762927575562172365337929,
  -23.48641564290853341361596821616234280392,
  1.362322948908907509911149920532561575254,
  -.7051467367205771043993968232310964220061e-2,
  0,
  0
};

typedef struct
{
  double *k[10];
  double *ytmp;
  double *y0;
}
vern7_state_t;

static void *
vern7_alloc (size_t dim)
{
  vern7_state_t *state = (vern7_state_t *) malloc (sizeof (vern7_state_t));
  int i, j;

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for vern7_state", GSL_ENOMEM);
    }

  state->ytmp = (double *) malloc (dim * sizeof (double));

  if (state->ytmp == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for ytmp", GSL_ENOMEM);
    }

  state->y0 = (double *) malloc (dim * sizeof (double));

  if (state->y0 == 0)
    {
      free (state->ytmp);
      free (state);
      GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
    }

  for (i = 0; i < 10; i++)
    {
      state->k[i] = (double *) malloc (dim * sizeof (double));

      if (state->k[i] == 0)
        {
          for (j = 0; j < i; j++)
            {
              free (state->k[j]);
            }
          free (state->y0);
          free (state->ytmp);
          free (state);
          GSL_ERROR_NULL ("failed to allocate space for k's", GSL_ENOMEM);
        }
    }

  return state;
}


static int
vern7_apply (void *vstate,
             size_t dim,
             double t,
             double h,
             double y[],
             double yerr[],
             const double dydt_in[],
             double dydt_out[], const gsl_odeiv2_system *sys)
{
  vern7_state_t *state = (vern7_state_t *) vstate;

  size_t i;

  double *const ytmp = state->ytmp;
  double *const y0 = state->y0;
  /* Note that k1 is stored in state->k[0] due to zero-based indexing */
  double *const k1 = state->k[0];
  double *const k2 = state->k[1];
  double *const k3 = state->k[2];
  double *const k4 = state->k[3];
  double *const k5 = state->k[4];
  double *const k6 = state->k[5];
  double *const k7 = state->k[6];
  double *const k8 = state->k[7];
  double *const k9 = state->k[8];
  double *const k10 = state->k[9];

  DBL_MEMCPY (y0, y, dim);

  /* k1 step */
  if (dydt_in != NULL)
    {
      DBL_MEMCPY (k1, dydt_in, dim);
    }
  else
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t, y, k1);

      if (s != GSL_SUCCESS)
        {
          return s;
        }
    }

  /* k2 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + b21 * h * k1[i];
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[0] * h, ytmp, k2);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k3 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b3[0] * k1[i] + b3[1] * k2[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[1] * h, ytmp, k3);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k4 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b4[0] * k1[i] + b4[2] * k3[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[2] * h, ytmp, k4);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k5 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b5[0] * k1[i] + b5[2] * k3[i] + b5[3] * k4[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[3] * h, ytmp, k5);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k6 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b6[0] * k1[i] + b6[2] * k3[i] + b6[3] * k4[i] +
                    b6[4] * k5[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[4] * h, ytmp, k6);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k7 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b7[0] * k1[i] + b7[2] * k3[i] + b7[3] * k4[i] +
                    b7[4] * k5[i] + b7[5] * k6[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[5] * h, ytmp, k7);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k8 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b8[0] * k1[i] + b8[2] * k3[i] + b8[3] * k4[i] +
                    b8[4] * k5[i] + b8[5] * k6[i] + b8[6] * k7[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[6] * h, ytmp, k8);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k9 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b9[0] * k1[i] + b9[2] * k3[i] + b9[3] * k4[i] +
                    b9[4] * k5[i] + b9[5] * k6[i] + b9[6] * k7[i] +
                    b9[7] * k8[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k9);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k10 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b10[0] * k1[i] + b10[2] * k3[i] + b10[3] * k4[i] +
                    b10[4] * k5[i] + b10[5] * k6[i] + b10[6] * k7[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k10);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* final sum  */
  for (i = 0; i < dim; i++)
    {
      const double ksum7 =
        Abar[0] * k1[i] + Abar[3] * k4[i] + Abar[4] * k5[i] +
        Abar[5] * k6[i] + Abar[6] * k7[i] + Abar[7] * k8[i] + Abar[8] * k9[i];
      y[i] += h * ksum7;
    }

  /* Evaluate dydt_out[]. */

  if (dydt_out != NULL)
    {
      int s = GSL_ODEIV_FN_EVAL (sys, t + h, y, dydt_out);

      if (s != GSL_SUCCESS)
        {
          /* Restore initial values */
          DBL_MEMCPY (y, y0, dim);
          return s;
        }
    }

  /* error estimate */
  for (i = 0; i < dim; i++)
    {
      const double ksum7 =
        Abar[0] * k1[i] + Abar[3] * k4[i] + Abar[4] * k5[i] +
        Abar[5] * k6[i] + Abar[6] * k7[i] + Abar[7] * k8[i] + Abar[8] * k9[i];
      const double ksum6 =
        A[0] * k1[i] + A[3] * k4[i] + A[4] * k5[i] + A[5] * k6[i] +
        A[6] * k7[i] + A[9] * k10[i];
      yerr[i] = h * (ksum6 - ksum7);
    }

  return GSL_SUCCESS;
}

static int
vern7_reset (void *vstate, size_t dim)
{
  vern7_state_t *state = (vern7_state_t *) vstate;

  int i;

  for (i = 0; i < 10; i++)
    {
      DBL_ZERO_MEMSET (state->k[i], dim);
    }

  DBL_ZERO_MEMSET (state->y0, dim);
  DBL_ZERO_MEMSET (state->ytmp, dim);

  return GSL_SUCCESS;
}

static unsigned int
vern7_order (void *vstate)
{
  vern7_state_t *state = (vern7_state_t *) vstate;
  state = 0;                    /* prevent warnings about unused parameters */
  return 7;
}

static void
vern7_free (void *vstate)
{
  vern7_state_t *state = (vern7_state_t *) vstate;
  int i;

  for (i = 0; i < 10; i++)
    {
      free (state->k[i]);
    }
  free (state->y0);
  free (state->ytmp);
  free (state);
}

static const gsl_odeiv2_step_type vern7_type = { "vern7",       /* name */
  1,                            /* can use dydt_in */
  1,                            /* gives exact dydt_out */
  &vern7_alloc,
  &vern7_apply,
  &stepper_set_driver_null,
  &vern7_reset,
  &vern7_order,
  &vern7_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_vern7 = &vern7_type;
