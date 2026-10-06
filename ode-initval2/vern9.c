/* ode-initval2/vern9.c
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

/* Runge-Kutta 9(8), Verner
 *
 * J.H. Verner, Numerically optimal Runge--Kutta pairs 
 * with interpolants. Numerical Algorithms, 53, (2010) pp. 383--396
 *
 * https://www.sfu.ca/~jverner/
 * "slightly more efficient" from April 2024.
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "odeiv_util.h"
#include "step_utils.c"

/* Vern9 constants */

static const double Abar[] = {
  .1500669014979724795766288712377040980022e-1,
  0,
  0,
  0,
  0,
  0,
  0,
  -1.055180992746381278594381685080184474710,
  .2384947263782183112638140557774851953952,
  .1288151774282991354622515887144081172371,
  .2276623111046215614614917896665746508574,
  1.229532587437517443160321815173576421896,
  .4624976662810383487397308275256128637703e-1,
  .1386196319366293903628308621201843492468,
  .3080010168319435405203560375162404389992e-1,
  0
};

static const double A[] = {
  .1897210532481101330735875918987801428423e-1,
  0,
  0,
  0,
  0,
  0,
  0,
  3.408110314549493848404398964228060776034,
  .1260323883820920906560270507988839661845,
  .1188375063451149770930378540709095182997,
  .2491041997838687569190073177760326138188,
  -3.269966219928978218713853055139116510335,
  .3023798100228882907409723963501699715077,
  0,
  0,
  .4652989552070924159305071272518165020637e-1
};

static const double ah[] = {
  .3571e-1,
  .9906028091267414072062290548990445721420e-1,
  .1485904213690112110809343582348566858213,
  .6134,
  .2327359473605626756631680289030288192566,
  .5538640526394373243368319710969711807434,
  .6555,
  .491625,
  .6858e-1,
  .253,
  .6620641795412045944786226689660879617864,
  .8309,
  .8998
};

static const double b21 = .3571e-1;

static const double b3[] = {
  -.3833735636677017025757228807792426014327e-1,
  .1373976372794443109781951935678287173575
};

static const double b4[] = {
  .3714760534225280277023358955871417145532e-1,
  0,
  .1114428160267584083107007686761425143660
};

static const double b5[] = {
  2.674764429871505119194043347886848033362,
  0,
  -9.982382134885293836441318015951354323506,
  7.921017705013788717247274668064506290145
};

static const double b6[] = {
  .5242104050577351069841401193318449756747e-1,
  0,
  0,
  .1796911189175953081279466649977364298948,
  .6237879371938568368073519721078917943307e-3
};

static const double b7[] = {
  .1592492223647632083060991150129732726703,
  0,
  0,
  -.4298429877241087508189185066909803751595,
  .6665266542726088051243445942394404211672e-1,
  .7578051525715219863372169033510342411158
};

static const double b8[] = {
  .7283333333333333333333333333333333333333e-1,
  0,
  0,
  0,
  0,
  .3359344590665103678713422141936031057620,
  .2467322076001562987953244524730635609046
};

static const double b9[] = {
  .7297558593750000000000000000000000000000e-1,
  0,
  0,
  0,
  0,
  .3348009729699333533129043555243262838993,
  .1184158239050666466870956444756737161007,
  -.3456738281250000000000000000000000000000e-1
};

static const double b10[] = {
  .4911213663452096382929799914888692571185e-1,
  0,
  0,
  0,
  0,
  .3983857361308652347066808423646954244674e-1,
  .1069675288939354812077139476183615284408,
  -.2174259165458647598455650055976622663086e-1,
  -.1055956474869564925231235304439517699685
};

static const double b11[] = {
  -.2707988818641280489175862864774930636734e-1,
  0,
  0,
  0,
  0,
  .3330000000000000000000000000000000000000e-1,
  -.1645526070036057075627040519129748090051,
  .3428266306497389946775103799535438271025e-1,
  .1585264064439221046855989221358760453749,
  .2185234256811225083011127204294936872873
};

static const double b12[] = {
  .5584657769108862727526281006220448941548e-1,
  0,
  0,
  0,
  0,
  .9166533166672539087479545790553151183505e-1,
  .2392399655523627049569179599100121305027,
  .1023834712248414879740316672186561113729e-1,
  -.2679331322859542442991625976085644567159e-2,
  .4235624181474284646680744174197106205439e-1,
  .2253970470166604185504274586005888014087
};

static const double b13[] = {
  -.4802510512725195893103945790000307739284,
  0,
  0,
  0,
  0,
  -6.359610162555930097843606206238737368598,
  -.2762313898040841385163666588363199062208,
  -6.500796633979846747011407331804389650124,
  .5734765877040956875980719881692630407295,
  1.347125994868138849716719646646049357523,
  5.936840409706221306933249144503479736029,
  6.590346245333924728433733996560685564589
};

static const double b14[] = {
  .3307533067671401079975082476263855100019,
  0,
  0,
  0,
  0,
  5.956207776829962111958631943398154563177,
  -.4868316400481527952968983900652017328476,
  4.462055288206771191616083778225319083232,
  .7410258231442071778525406966672954745201,
  -.7118192034575913119318346652024763899638,
  -5.454619594516665440649148499095762263275,
  -4.140803729244709693305868672684077643030,
  .2038319723190386517589855611303633981843
};

static const double b15[] = {
  -.5847111122998944389651719860007729221315,
  0,
  0,
  0,
  0,
  -12.41268417116267068926060969930886386670,
  1.360245445660928146977003413153919520189,
  -22.42610531111868229005965344123311160860,
  -.8828857055865458252669099251378431470909,
  1.770155128538230466813900572142546329490,
  12.15809651918533877218739262890752294243,
  22.23037520407760703834156132941550000589,
  -.6634483760201249006317401316871863808959,
  .4509623787258137198642272397482891274306
};

static const double b16[] = {
  1.940575549810648717395482898478643635949,
  0,
  0,
  0,
  0,
  21.97798408114556310164832153356570429956,
  .8230747326984728145128172720328075817426,
  68.16441683626354817206843139989369741474,
  -3.117097463620266666801525529574905518624,
  -4.568841021822439620303730003755793384193,
  -18.74190987126264955904919353670543412590,
  -66.57711839637831878350036856645874492472,
  1.098915553165441824029764532524025021446,
  0,
  0
};

typedef struct
{
  double *k[16];
  double *ytmp;
  double *y0;
}
vern9_state_t;

static void *
vern9_alloc (size_t dim)
{
  vern9_state_t *state = (vern9_state_t *) malloc (sizeof (vern9_state_t));
  int i, j;

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for vern9_state", GSL_ENOMEM);
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

  for (i = 0; i < 16; i++)
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
vern9_apply (void *vstate,
             size_t dim,
             double t,
             double h,
             double y[],
             double yerr[],
             const double dydt_in[],
             double dydt_out[], const gsl_odeiv2_system *sys)
{
  vern9_state_t *state = (vern9_state_t *) vstate;

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
  double *const k14 = state->k[13];
  double *const k15 = state->k[14];
  double *const k16 = state->k[15];

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
      ytmp[i] = y[i] + h * (b8[0] * k1[i] + b8[5] * k6[i] + b8[6] * k7[i]);
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
        y[i] + h * (b9[0] * k1[i] + b9[5] * k6[i] + b9[6] * k7[i] +
                    b9[7] * k8[i]);
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
        y[i] + h * (b10[0] * k1[i] + b10[5] * k6[i] + b10[6] * k7[i] +
                    b10[7] * k8[i] + b10[8] * k9[i]);
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
        y[i] + h * (b11[0] * k1[i] + b11[5] * k6[i] + b11[6] * k7[i] +
                    b11[7] * k8[i] + b11[8] * k9[i] + b11[9] * k10[i]);
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
        y[i] + h * (b12[0] * k1[i] + b12[5] * k6[i] + b12[6] * k7[i] +
                    b12[7] * k8[i] + b12[8] * k9[i] + b12[9] * k10[i] +
                    b12[10] * k11[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[10] * h, ytmp, k12);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k13 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b13[0] * k1[i] + b13[5] * k6[i] + b13[6] * k7[i] +
                    b13[7] * k8[i] + b13[8] * k9[i] + b13[9] * k10[i] +
                    b13[10] * k11[i] + b13[11] * k12[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[11] * h, ytmp, k13);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k14 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b14[0] * k1[i] + b14[5] * k6[i] + b14[6] * k7[i] +
                    b14[7] * k8[i] + b14[8] * k9[i] + b14[9] * k10[i] +
                    b14[10] * k11[i] + b14[11] * k12[i] + b14[12] * k13[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + ah[12] * h, ytmp, k14);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k15 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b15[0] * k1[i] + b15[5] * k6[i] + b15[6] * k7[i] +
                    b15[7] * k8[i] + b15[8] * k9[i] + b15[9] * k10[i] +
                    b15[10] * k11[i] + b15[11] * k12[i] + b15[12] * k13[i] +
                    b15[13] * k14[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k15);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k16 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b16[0] * k1[i] + b16[5] * k6[i] + b16[6] * k7[i] +
                    b16[7] * k8[i] + b16[8] * k9[i] + b16[9] * k10[i] +
                    b16[10] * k11[i] + b16[11] * k12[i] + b16[12] * k13[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k16);

    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* final sum  */
  for (i = 0; i < dim; i++)
    {
      const double ksum9 =
        Abar[0] * k1[i] + Abar[7] * k8[i] + Abar[8] * k9[i] +
        Abar[9] * k10[i] + Abar[10] * k11[i] + Abar[11] * k12[i] +
        Abar[12] * k13[i] + Abar[13] * k14[i] + Abar[14] * k15[i];
      y[i] += h * ksum9;
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
      const double ksum9 =
        Abar[0] * k1[i] + Abar[7] * k8[i] + Abar[8] * k9[i] +
        Abar[9] * k10[i] + Abar[10] * k11[i] + Abar[11] * k12[i] +
        Abar[12] * k13[i] + Abar[13] * k14[i] + Abar[14] * k15[i];
      const double ksum8 =
        A[0] * k1[i] + A[7] * k8[i] + A[8] * k9[i] + A[9] * k10[i] +
        A[10] * k11[i] + A[11] * k12[i] + A[12] * k13[i] + A[15] * k16[i];
      yerr[i] = h * (ksum8 - ksum9);
    }

  return GSL_SUCCESS;
}

static int
vern9_reset (void *vstate, size_t dim)
{
  vern9_state_t *state = (vern9_state_t *) vstate;

  int i;

  for (i = 0; i < 16; i++)
    {
      DBL_ZERO_MEMSET (state->k[i], dim);
    }

  DBL_ZERO_MEMSET (state->y0, dim);
  DBL_ZERO_MEMSET (state->ytmp, dim);

  return GSL_SUCCESS;
}

static unsigned int
vern9_order (void *vstate)
{
  vern9_state_t *state = (vern9_state_t *) vstate;
  state = 0;                    /* prevent warnings about unused parameters */
  return 9;
}

static void
vern9_free (void *vstate)
{
  vern9_state_t *state = (vern9_state_t *) vstate;
  int i;

  for (i = 0; i < 16; i++)
    {
      free (state->k[i]);
    }
  free (state->y0);
  free (state->ytmp);
  free (state);
}

static const gsl_odeiv2_step_type vern9_type = { "vern9",       /* name */
  1,                            /* can use dydt_in */
  1,                            /* gives exact dydt_out */
  &vern9_alloc,
  &vern9_apply,
  &stepper_set_driver_null,
  &vern9_reset,
  &vern9_order,
  &vern9_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_vern9 = &vern9_type;
