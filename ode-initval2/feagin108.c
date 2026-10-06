/* ode-initval2/feagin108.c
 * 
 * Copyright (C) 2024 Christoph M. Schaefer 
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

/* Feagin Runge-Kutta 10(8) */

/* Reference: HIGH-ORDER EXPLICIT RUNGE-KUTTA
   METHODS USING m-SYMMETRY, T. Feagin, 
   Neural, Parallel, and Scientific Computations 20 (2012) 437-458 */

/* Author:  C.M. Schaefer
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "odeiv_util.h"
#include "step_utils.c"

/* Feagin RK10(8) constants, see https://sce.uhcl.edu/rungekutta/rk108.txt */

static const double a0 = 0.0;
static const double a1 =
  0.100000000000000000000000000000000000000000000000000000000000;
static const double a2 =
  0.539357840802981787532485197881302436857273449701009015505500;
static const double a3 =
  0.809036761204472681298727796821953655285910174551513523258250;
static const double a4 =
  0.309036761204472681298727796821953655285910174551513523258250;
static const double a5 =
  0.981074190219795268254879548310562080489056746118724882027805;
static const double a6 =
  0.833333333333333333333333333333333333333333333333333333333333;
static const double a7 =
  0.354017365856802376329264185948796742115824053807373968324184;
static const double a8 =
  0.882527661964732346425501486979669075182867844268052119663791;
static const double a9 =
  0.642615758240322548157075497020439535959501736363212695909875;
static const double a10 =
  0.357384241759677451842924502979560464040498263636787304090125;
static const double a11 =
  0.117472338035267653574498513020330924817132155731947880336209;
static const double a12 =
  0.833333333333333333333333333333333333333333333333333333333333;
static const double a13 =
  0.309036761204472681298727796821953655285910174551513523258250;
static const double a14 =
  0.539357840802981787532485197881302436857273449701009015505500;
static const double a15 =
  0.100000000000000000000000000000000000000000000000000000000000;
static const double a16 = 1.0;

static const double c0 =
  0.0333333333333333333333333333333333333333333333333333333333333;
static const double c1 =
  0.0250000000000000000000000000000000000000000000000000000000000;
static const double c2 =
  0.0333333333333333333333333333333333333333333333333333333333333;
static const double c4 =
  0.0500000000000000000000000000000000000000000000000000000000000;
static const double c6 =
  0.0400000000000000000000000000000000000000000000000000000000000;
static const double c8 =
  0.189237478148923490158306404106012326238162346948625830327194;
static const double c9 =
  0.277429188517743176508360262560654340428504319718040836339472;
static const double c10 =
  0.277429188517743176508360262560654340428504319718040836339472;
static const double c11 =
  0.189237478148923490158306404106012326238162346948625830327194;
static const double c12 =
  -0.0400000000000000000000000000000000000000000000000000000000000;
static const double c13 =
  -0.0500000000000000000000000000000000000000000000000000000000000;
static const double c14 =
  -0.0333333333333333333333333333333333333333333333333333333333333;
static const double c15 =
  -0.0250000000000000000000000000000000000000000000000000000000000;
static const double c16 =
  0.0333333333333333333333333333333333333333333333333333333333333;

static const double b10 =
  0.100000000000000000000000000000000000000000000000000000000000;
static const double b20 =
  -0.915176561375291440520015019275342154318951387664369720564660;
static const double b21 =
  1.45453440217827322805250021715664459117622483736537873607016;
static const double b30 =
  0.202259190301118170324681949205488413821477543637878380814562;
static const double b32 =
  0.606777570903354510974045847616465241464432630913635142443687;
static const double b40 =
  0.184024714708643575149100693471120664216774047979591417844635;
static const double b42 =
  0.197966831227192369068141770510388793370637287463360401555746;
static const double b43 =
  -0.0729547847313632629185146671595558023015011608914382961421311;
static const double b50 =
  0.0879007340206681337319777094132125475918886824944548534041378;
static const double b53 =
  0.410459702520260645318174895920453426088035325902848695210406;
static const double b54 =
  0.482713753678866489204726942976896106809132737721421333413261;
static const double b60 =
  0.0859700504902460302188480225945808401411132615636600222593880;
static const double b63 =
  0.330885963040722183948884057658753173648240154838402033448632;
static const double b64 =
  0.489662957309450192844507011135898201178015478433790097210790;
static const double b65 =
  -0.0731856375070850736789057580558988816340355615025188195854775;
static const double b70 =
  0.120930449125333720660378854927668953958938996999703678812621;
static const double b74 =
  0.260124675758295622809007617838335174368108756484693361887839;
static const double b75 =
  0.0325402621549091330158899334391231259332716675992700000776101;
static const double b76 =
  -0.0595780211817361001560122202563305121444953672762930724538856;
static const double b80 =
  0.110854379580391483508936171010218441909425780168656559807038;
static const double b85 =
  -0.0605761488255005587620924953655516875526344415354339234619466;
static const double b86 =
  0.321763705601778390100898799049878904081404368603077129251110;
static const double b87 =
  0.510485725608063031577759012285123416744672137031752354067590;
static const double b90 =
  0.112054414752879004829715002761802363003717611158172229329393;
static const double b95 =
  -0.144942775902865915672349828340980777181668499748506838876185;
static const double b96 =
  -0.333269719096256706589705211415746871709467423992115497968724;
static const double b97 =
  0.499269229556880061353316843969978567860276816592673201240332;
static const double b98 =
  0.509504608929686104236098690045386253986643232352989602185060;
static const double b100 =
  0.113976783964185986138004186736901163890724752541486831640341;
static const double b105 =
  -0.0768813364203356938586214289120895270821349023390922987406384;
static const double b106 =
  0.239527360324390649107711455271882373019741311201004119339563;
static const double b107 =
  0.397774662368094639047830462488952104564716416343454639902613;
static const double b108 =
  0.0107558956873607455550609147441477450257136782823280838547024;
static const double b109 =
  -0.327769124164018874147061087350233395378262992392394071906457;
static const double b110 =
  0.0798314528280196046351426864486400322758737630423413945356284;
static const double b115 =
  -0.0520329686800603076514949887612959068721311443881683526937298;
static const double b116 =
  -0.0576954146168548881732784355283433509066159287152968723021864;
static const double b117 =
  0.194781915712104164976306262147382871156142921354409364738090;
static const double b118 =
  0.145384923188325069727524825977071194859203467568236523866582;
static const double b119 =
  -0.0782942710351670777553986729725692447252077047239160551335016;
static const double b1110 =
  -0.114503299361098912184303164290554670970133218405658122674674;
static const double b120 =
  0.985115610164857280120041500306517278413646677314195559520529;
static const double b123 =
  0.330885963040722183948884057658753173648240154838402033448632;
static const double b124 =
  0.489662957309450192844507011135898201178015478433790097210790;
static const double b125 =
  -1.37896486574843567582112720930751902353904327148559471526397;
static const double b126 =
  -0.861164195027635666673916999665534573351026060987427093314412;
static const double b127 =
  5.78428813637537220022999785486578436006872789689499172601856;
static const double b128 =
  3.28807761985103566890460615937314805477268252903342356581925;
static const double b129 =
  -2.38633905093136384013422325215527866148401465975954104585807;
static const double b1210 =
  -3.25479342483643918654589367587788726747711504674780680269911;
static const double b1211 =
  -2.16343541686422982353954211300054820889678036420109999154887;
static const double b130 =
  0.895080295771632891049613132336585138148156279241561345991710;
static const double b132 =
  0.197966831227192369068141770510388793370637287463360401555746;
static const double b133 =
  -0.0729547847313632629185146671595558023015011608914382961421311;
static const double b135 =
  -0.851236239662007619739049371445966793289359722875702227166105;
static const double b136 =
  0.398320112318533301719718614174373643336480918103773904231856;
static const double b137 =
  3.63937263181035606029412920047090044132027387893977804176229;
static const double b138 =
  1.54822877039830322365301663075174564919981736348973496313065;
static const double b139 =
  -2.12221714704053716026062427460427261025318461146260124401561;
static const double b1310 =
  -1.58350398545326172713384349625753212757269188934434237975291;
static const double b1311 =
  -1.71561608285936264922031819751349098912615880827551992973034;
static const double b1312 =
  -0.0244036405750127452135415444412216875465593598370910566069132;
static const double b140 =
  -0.915176561375291440520015019275342154318951387664369720564660;
static const double b141 =
  1.45453440217827322805250021715664459117622483736537873607016;
static const double b144 =
  -0.777333643644968233538931228575302137803351053629547286334469;
static const double b146 =
  -0.0910895662155176069593203555807484200111889091770101799647985;
static const double b1412 =
  0.0910895662155176069593203555807484200111889091770101799647985;
static const double b1413 =
  0.777333643644968233538931228575302137803351053629547286334469;
static const double b150 =
  0.100000000000000000000000000000000000000000000000000000000000;
static const double b152 =
  -0.157178665799771163367058998273128921867183754126709419409654;
static const double b1514 =
  0.157178665799771163367058998273128921867183754126709419409654;
static const double b160 =
  0.181781300700095283888472062582262379650443831463199521664945;
static const double b161 =
  0.675000000000000000000000000000000000000000000000000000000000;
static const double b162 =
  0.342758159847189839942220553413850871742338734703958919937260;
static const double b164 =
  0.259111214548322744512977076191767379267783684543182428778156;
static const double b165 =
  -0.358278966717952089048961276721979397739750634673268802484271;
static const double b166 =
  -1.04594895940883306095050068756409905131588123172378489286080;
static const double b167 =
  0.930327845415626983292300564432428777137601651182965794680397;
static const double b168 =
  1.77950959431708102446142106794824453926275743243327790536000;
static const double b169 =
  0.100000000000000000000000000000000000000000000000000000000000;
static const double b1610 =
  -0.282547569539044081612477785222287276408489375976211189952877;
static const double b1611 =
  -0.159327350119972549169261984373485859278031542127551931461821;
static const double b1612 =
  -0.145515894647001510860991961081084111308650130578626404945571;
static const double b1613 =
  -0.259111214548322744512977076191767379267783684543182428778156;
static const double b1614 =
  -0.342758159847189839942220553413850871742338734703958919937260;
static const double b1615 =
  -0.675000000000000000000000000000000000000000000000000000000000;


typedef struct
{
  double *k[17];
  double *y0;
  double *ytmp;
}
feagin108_state_t;

static void *
feagin108_alloc (size_t dim)
{
  feagin108_state_t *state =
    (feagin108_state_t *) malloc (sizeof (feagin108_state_t));
  int i, j;

  if (state == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for feagin108_state",
                      GSL_ENOMEM);
    }

  state->ytmp = (double *) malloc (dim * sizeof (double));

  if (state->ytmp == 0)
    {
      free (state);
      GSL_ERROR_NULL ("failed to allocated space for ytmp", GSL_ENOMEM);
    }

  state->y0 = (double *) malloc (dim * sizeof (double));

  if (state->y0 == 0)
    {
      free (state->ytmp);
      free (state);
      GSL_ERROR_NULL ("failed to allocated space for y0", GSL_ENOMEM);
    }

  for (i = 0; i < 17; i++)
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
feagin108_apply (void *vstate,
                 size_t dim,
                 double t,
                 double h,
                 double y[],
                 double yerr[],
                 const double dydt_in[],
                 double dydt_out[], const gsl_odeiv2_system *sys)
{
  feagin108_state_t *state = (feagin108_state_t *) vstate;
  size_t i;

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
  double *const k17 = state->k[16];
  double *const y0 = state->y0;
  double *const ytmp = state->ytmp;


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
      ytmp[i] = y[i] + h * b10 * k1[i];
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a1 * h, ytmp, k2);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k3 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b20 * k1[i] + b21 * k2[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a2 * h, ytmp, k3);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k4 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b30 * k1[i] + b32 * k3[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a3 * h, ytmp, k4);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k5 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b40 * k1[i] + b42 * k3[i] + b43 * k4[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a4 * h, ytmp, k5);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k6 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b50 * k1[i] + b53 * k4[i] + b54 * k5[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a5 * h, ytmp, k6);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k7 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b60 * k1[i] + b63 * k4[i] + b64 * k5[i] + b65 * k6[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a6 * h, ytmp, k7);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k8 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b70 * k1[i] + b74 * k5[i] + b75 * k6[i] + b76 * k7[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a7 * h, ytmp, k8);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k9 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b80 * k1[i] + b85 * k6[i] + b86 * k7[i] + b87 * k8[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a8 * h, ytmp, k9);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k10 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b90 * k1[i] + b95 * k6[i] + b96 * k7[i] + b97 * k8[i] +
                    b98 * k9[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a9 * h, ytmp, k10);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k11 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b100 * k1[i] + b105 * k6[i] + b106 * k7[i] +
                    b107 * k8[i] + b108 * k9[i] + b109 * k10[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a10 * h, ytmp, k11);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k12 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b110 * k1[i] + b115 * k6[i] + b116 * k7[i] +
                    b117 * k8[i] + b118 * k9[i] + b119 * k10[i] +
                    b1110 * k11[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a11 * h, ytmp, k12);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k13 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b120 * k1[i] + b123 * k4[i] + b124 * k5[i] +
                    b125 * k6[i] + b126 * k7[i] + b127 * k8[i] +
                    b128 * k9[i] + b129 * k10[i] + b1210 * k11[i] +
                    b1211 * k12[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a12 * h, ytmp, k13);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k14 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b130 * k1[i] + b132 * k3[i] + b133 * k4[i] +
                    b135 * k6[i] + b136 * k7[i] + b137 * k8[i] +
                    b138 * k9[i] + b139 * k10[i] + b1310 * k11[i] +
                    b1311 * k12[i] + b1312 * k13[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a13 * h, ytmp, k14);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k15 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b140 * k1[i] + b141 * k2[i] + b144 * k5[i] +
                    b146 * k7[i] + b1412 * k13[i] + b1413 * k14[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a14 * h, ytmp, k15);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k16 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] = y[i] + h * (b150 * k1[i] + b152 * k3[i] + b1514 * k15[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a15 * h, ytmp, k16);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* k17 step */
  for (i = 0; i < dim; i++)
    {
      ytmp[i] =
        y[i] + h * (b160 * k1[i] + b161 * k2[i] + b162 * k3[i] +
                    b164 * k5[i] + b165 * k6[i] + b166 * k7[i] +
                    b167 * k8[i] + b168 * k9[i] + b169 * k10[i] +
                    b1610 * k11[i] + b1611 * k12[i] + b1612 * k13[i] +
                    b1613 * k14[i] + b1614 * k15[i] + b1615 * k16[i]);
    }

  {
    int s = GSL_ODEIV_FN_EVAL (sys, t + a16 * h, ytmp, k17);
    if (s != GSL_SUCCESS)
      {
        return s;
      }
  }

  /* final sum */
  for (i = 0; i < dim; i++)
    {
      const double d_i = c0 * k1[i] + c1 * k2[i] + c2 * k3[i] + c4 * k5[i]
        + c6 * k7[i] + c8 * k9[i]
        + c9 * k10[i] + c10 * k11[i] + c11 * k12[i] + c12 * k13[i]
        + c13 * k14[i] + c14 * k15[i] + c15 * k16[i] + c16 * k17[i];

      y[i] += h * d_i;
    }


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


  /* difference between 12th and 14th order */
  for (i = 0; i < dim; i++)
    {
      yerr[i] = 1. / 360. * (k2[i] - k16[i]) * h;
    }

  return GSL_SUCCESS;

}


static int
feagin108_reset (void *vstate, size_t dim)
{
  feagin108_state_t *state = (feagin108_state_t *) vstate;

  int i;

  for (i = 0; i < 17; i++)
    {
      DBL_ZERO_MEMSET (state->k[i], dim);
    }

  DBL_ZERO_MEMSET (state->ytmp, dim);
  DBL_ZERO_MEMSET (state->y0, dim);

  return GSL_SUCCESS;
}

static unsigned int
feagin108_order (void *vstate)
{
  feagin108_state_t *state = (feagin108_state_t *) vstate;
  state = 0;                    /* prevent warnings about unused parameters */
  return 10;
}

static void
feagin108_free (void *vstate)
{
  feagin108_state_t *state = (feagin108_state_t *) vstate;
  int i;

  for (i = 0; i < 17; i++)
    {
      free (state->k[i]);
    }

  free (state->ytmp);
  free (state->y0);
  free (state);
}

static const gsl_odeiv2_step_type feagin108_type = { "feagin108",       /* name */
  1,                            /* can use dydt_in */
  1,                            /* gives exact dydt_out */
  &feagin108_alloc,
  &feagin108_apply,
  &stepper_set_driver_null,
  &feagin108_reset,
  &feagin108_order,
  &feagin108_free
};

const gsl_odeiv2_step_type *gsl_odeiv2_step_feagin108 = &feagin108_type;
