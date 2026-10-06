/* ode-initval2/vern8.c
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

/* Runge-Kutta 8(7), Verner
 *
 * J.H. Verner, Numerically optimal Runge--Kutta pairs 
 * with interpolants. Numerical Algorithms, 53, (2010) pp. 383--396
 *
 * https://www.sfu.ca/~jverner/
 * "considerably more efficient" from May 2024.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "odeiv_util.h"
#include "step_utils.c"

/* Vern8 constants */

static const double ah[] = {
  .92662e-1,
  .1312230361754017604780747799402406525075,
  .1968345542631026407171121699103609787612,
  .427173,
  .485972,
  .161915,
  .985468,
  .9626977348604540392115664231947926180050,
  .99626,
  .997947
};

static const double b21 = .92662e-1;

static const double b3[] = {
  .3830746548250284242039554953085778876548e-1,
  .9291557069289891805767923040938286374199e-1
};

static const double b4[] = {
  .4920863856577566017927804247759024469030e-1,
  0,
  .1476259156973269805378341274327707340709
};

static const double b5[] = {
  .2743076085702486894953667699331084997155,
  0,
  -.9319887203102656329703895298805780048911,
  1.084854111740016943475022759947469505176
};

static const double b6[] = {
  .6461852970939692178473977038055776659235e-1,
  0,
  0,
  .2687629213368923356417418162215129378455,
  .1525905489537107425735184133979292955622
};

static const double b7[] = {
  .7189155819773216802747111989237192757241e-1,
  0,
  0,
  .1221265783362549661095644262553368220549,
  -.7943550859198561207449556225926957458888e-1,
  .4733237205799847793746001611156082496161e-1
};

static const double b8[] = {
  -6.073603893714328779581044654669389532243,
  0,
  0,
  -73.8956,
  11.93985370695273926305714852083567777565,
  -3.839251541405054537968455811084326767757,
  72.85406972816664405449235194491803852435
};

static const double b9[] = {
  -4.868640079323569115532909137641679122100,
  0,
  0,
  -59.18572799975646020086112235306397792921,
  9.230819319232425236436363261556755728972,
  -2.676847914962525780976208768457285096420,
  58.45720009994685754009028962323039746457,
  .5894309723726360055153797570581572198272e-2
};

static const double b10[] = {
  -6.689861899320853351893486767784362717327,
  0,
  0,
  -81.44271004053111286646869388606624375878,
  13.36778825698397107436768361456719138955,
  -4.470777638416181156550300441220009175768,
  80.23321392161410397716314034458127569922,
  -.1313638336212181564579824084323466067556e-1,
  .1174378303219413902745537676538322377075e-1
};

static const double b11[] = {
  -6.788841955800464091202279324616261399375,
  0,
  0,
  -82.65639855934828888595908426400236814314,
  13.59973921874899036529372991555792312537,
  -4.574464055350503753023585975406886448931,
  81.41943207216075927716035509139094411649,
  -.1416248014826418017596237345824563213708e-1,
  .1375441580835227405704608471667945126391e-1,
  -.1111656070581006150219154181785069541859e-2
};

static const double b12[] = {
  -6.910189846402485729483801895952800612193,
  0,
  0,
  -84.14495154176748682468399818934446849309,
  13.88512122378983816888937012217322612857,
  -4.702458788144493296650978281402083070191,
  82.87411451529241610017736751027256505872,
  -.1645498337198780129448623133684652608689e-1,
  .1644663972162521365153657249183964172873e-1,
  .4275449370796530995842729680318146579263e-2,
  -.5902668488222361600852336581750274034670e-2
};

static const double b13[] = {
  -6.911973921198979615960353484856419992328,
  0,
  0,
  -84.16635595878781036984379688648536865215,
  13.88834627565582007275122755367241993143,
  -4.703463178409702575934171709526131267400,
  82.89518622207404915885492275044886922215,
  -.1020345016228260287805622619295788864042e-1,
  .1427900423230391471526298444108039935137e-1,
  -.5814993403397981705034981501491752405078e-2
};

static const double Abar[] = {
  .4625543159712467285354070519930680076661e-1,
  0,
  0,
  0,
  0,
  .3706666165521011182439275381303388440188,
  .2590440824552746577195309846039127860158,
  -679.9841468175039046601229652340421033215,
  49.89161129042053159104301060910837813887,
  10271.23522213731241388784467688648118861,
  -14782.19660635689728059010570422698640618,
  5141.377953616063739322523982737505384322,
  0
};

static const double A[] = {
  .4638504234365210644214797353760063769606e-1,
  0,
  0,
  0,
  0,
  .3725767681581196020016675337652919509451,
  .2585685495121687311414935155287094016045,
  -147.4950767589265301611858882107047892762,
  23.84362712644587506762535723842314335392,
  347.4264166730550802334952839143352280200,
  0,
  0,
  -223.4524974005883655795200619648851840879
};

typedef struct
{
  double *k[13];
  double *ytmp;
  double *y0;
}
vern8_state_t;

static void *
vern8_alloc (size_t dim)
{
  vern8_state_t *state = (vern8_state_t *) malloc (sizeof (vern8_state_t));
  int i, j;

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for vern8_state", GSL_ENOMEM);
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

  for (i = 0; i < 13; i++)
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
vern8_apply (void *vstate,
             size_t dim,
             double t,
             double h,
             double y[],
             double yerr[],
             const double dydt_in[],
             double dydt_out[], const gsl_odeiv2_system *sys)
{
  vern8_state_t *state = (vern8_state_t *) vstate;

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
  double *const k11 = state->k[10];
  double *const k12 = state->k[11];
  double *const k13 = state->k[12];

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
      ytmp[i] = y[i] + h * (b6[0] * k1[i] + b6[3] * k4[i] + b6[4] * k5[i]);
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
        y[i] + h * (b7[0] * k1[i] + b7[3] * k4[i] + b7[4] * k5[i] +
                    b7[5] * k6[i]);
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
        y[i] + h * (b8[0] * k1[i] + b8[3] * k4[i] + b8[4] * k5[i] +
                    b8[5] * k6[i] + b8[6] * k7[i]);
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
        y[i] + h * (b9[0] * k1[i] + b9[3] * k4[i] + b9[4] * k5[i] +
                    b9[5] * k6[i] + b9[6] * k7[i] + b9[7] * k8[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[7] * h, ytmp, k9);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k10 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b10[0] * k1[i] + b10[3] * k4[i] + b10[4] * k5[i] +
                    b10[5] * k6[i] + b10[6] * k7[i] + b10[7] * k8[i] +
                    b10[8] * k9[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[8] * h, ytmp, k10);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k11 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b11[0] * k1[i] + b11[3] * k4[i] + b11[4] * k5[i] +
                    b11[5] * k6[i] + b11[6] * k7[i] + b11[7] * k8[i] +
                    b11[8] * k9[i] + b11[9] * k10[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[9] * h, ytmp, k11);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k12 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b12[0] * k1[i] + b12[3] * k4[i] + b12[4] * k5[i] +
                    b12[5] * k6[i] + b12[6] * k7[i] + b12[7] * k8[i] +
                    b12[8] * k9[i] + b12[9] * k10[i] + b12[10] * k11[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k12);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k13 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b13[0] * k1[i] + b13[3] * k4[i] + b13[4] * k5[i] +
                    b13[5] * k6[i] + b13[6] * k7[i] + b13[7] * k8[i] +
                    b13[8] * k9[i] + b13[9] * k10[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k13);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* final sum  */
  for (i = 0; i < dim; i++)
    {
      const double ksum8 =
        Abar[0] * k1[i] + Abar[5] * k6[i] + Abar[6] * k7[i] +
        Abar[7] * k8[i] + Abar[8] * k9[i] + Abar[9] * k10[i] +
        Abar[10] * k11[i] + Abar[11] * k12[i];
      y[i] += h * ksum8;
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
      const double ksum8 =
        Abar[0] * k1[i] + Abar[5] * k6[i] + Abar[6] * k7[i] +
        Abar[7] * k8[i] + Abar[8] * k9[i] + Abar[9] * k10[i] +
        Abar[10] * k11[i] + Abar[11] * k12[i];
      const double ksum7 =
        A[0] * k1[i] + A[5] * k6[i] + A[6] * k7[i] + A[7] * k8[i] +
        A[8] * k9[i] + A[9] * k10[i] + A[12] * k13[i];
      yerr[i] = h * (ksum7 - ksum8);
    }

  return GSL_SUCCESS;
}

static int
vern8_reset (void *vstate, size_t dim)
{
  vern8_state_t *state = (vern8_state_t *) vstate;

  int i;

  for (i = 0; i < 13; i++)
    {
      DBL_ZERO_MEMSET (state->k[i], dim);
    }

  DBL_ZERO_MEMSET (state->y0, dim);
  DBL_ZERO_MEMSET (state->ytmp, dim);

  return GSL_SUCCESS;
}

static unsigned int
vern8_order (void *vstate)
{
  vern8_state_t *state = (vern8_state_t *) vstate;
  state = 0;                    /* prevent warnings about unused parameters */
  return 8;
}

static void
vern8_free (void *vstate)
{
  vern8_state_t *state = (vern8_state_t *) vstate;
  int i;

  for (i = 0; i < 13; i++)
    {
      free (state->k[i]);
    }
  free (state->y0);
  free (state->ytmp);
  free (state);
}

static const gsl_odeiv2_step_type vern8_type = { "vern8",       /* name */
  1,                            /* can use dydt_in */
  1,                            /* gives exact dydt_out */
  &vern8_alloc,
  &vern8_apply,
  &stepper_set_driver_null,
  &vern8_reset,
  &vern8_order,
  &vern8_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_vern8 = &vern8_type;
