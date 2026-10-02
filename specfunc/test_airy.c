/* specfunc/test_airy.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"


int test_airy(void)
{
  int s = 0;
  int m = GSL_MODE_DEFAULT;
  gsl_sf_result r;

  /** functions */

  TEST_SF(s, gsl_sf_airy_Ai_e, (-500.0, m, &r),              0.0725901201040411396, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-5.0, m, &r),                0.3507610090241142,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-0.3000000000000094, m, &r), 0.4309030952855831,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (0.6999999999999907, m, &r),  0.1891624003981519,    TEST_TOL0, GSL_SUCCESS);

  /*  This original value seemed to be slightly inaccurate in the last place.
      I recomputed it with pari to get the new value which end in 885 
      instead of 882 */
  /*
    TEST_SF(s, gsl_sf_airy_Ai_e, (1.649999999999991, m, &r),   0.05831058618720882,   TEST_TOL0, GSL_SUCCESS);
    */
  
  TEST_SF(s, gsl_sf_airy_Ai_e, (1.649999999999991, m, &r),   0.0583105861872088521,   TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Ai_e, (2.54999999999999, m, &r),    0.01446149513295428,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (3.499999999999987, m, &r),   0.002584098786989702,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (5.39999999999998, m, &r),    4.272986169411866e-05, TEST_TOL0, GSL_SUCCESS);
  
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (-5.0, m, &r),                  0.3507610090241142, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (0.6999999999999907, m, &r), 0.2795125667681217, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (1.649999999999991, m, &r),  0.2395493001442741, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (2.54999999999999, m, &r),   0.2183658595899388, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (3.499999999999987, m, &r),  0.2032920808163519, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_scaled_e, (5.39999999999998, m, &r),   0.1836050093282229, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Bi_e, (-500.0, m, &r),             -0.094688570132991028, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (-5.0, m, &r),               -0.1383691349016005,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (0.6999999999999907, m, &r),  0.9733286558781599,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (1.649999999999991, m, &r),   2.196407956850028,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (2.54999999999999, m, &r),    6.973628612493443,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (3.499999999999987, m, &r),   33.05550675461069,    TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_e, (5.39999999999998, m, &r),    1604.476078241272,    TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (-5.0, m, &r),                  -0.1383691349016005, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (0.6999999999999907, m, &r),  0.6587080754582302, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (1.649999999999991, m, &r),   0.5346449995597539, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (2.54999999999999, m, &r),    0.461835455542297,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (3.499999999999987, m, &r),   0.4201771882353061, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_scaled_e, (5.39999999999998, m, &r),    0.3734050675720473, TEST_TOL0, GSL_SUCCESS);


  /** derivatives */

  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (-5.0, m, &r),                 0.3271928185544435,       TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (-0.5500000000000094, m, &r), -0.1914604987143629,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (0.4999999999999906, m, &r),  -0.2249105326646850,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (1.899999999999992, m, &r),   -0.06043678178575718,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (3.249999999999988, m, &r),   -0.007792687926790889,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (5.199999999999981, m, &r),   -0.0001589434526459543, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (-5.0, m, &r),                0.3271928185544435, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (0.5499999999999906, m, &r), -0.2874057279170166, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (1.499999999999991, m, &r),  -0.3314199796863637, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (2.49999999999999, m, &r),   -0.3661089384751620, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (3.649999999999986, m, &r),  -0.3974033831453963, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_scaled_e, (6.299999999999977, m, &r),  -0.4508799189585947, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (-5.0, m, &r),                0.778411773001899,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (-0.5500000000000094, m, &r), 0.5155785358765014, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (0.4999999999999906, m, &r),  0.5445725641405883, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (1.899999999999992, m, &r),   3.495165862891568,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (3.249999999999988, m, &r),   36.55485149250338,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (5.199999999999981, m, &r),   2279.748293583233,  TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (-5.0, m, &r),               0.778411773001899,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (0.5499999999999906, m, &r), 0.4322811281817566, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (1.499999999999991, m, &r),  0.5542307563918037, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (2.49999999999999, m, &r),   0.6755384441644985, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (3.649999999999986, m, &r),  0.7613959373000228, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_scaled_e, (6.299999999999977, m, &r),  0.8852064139737571, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (2, &r),  -4.087949444130970617, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (50,   &r), -38.02100867725525443, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (100, &r),  -60.45555727411669871, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (110, &r),  -64.43135670991324811, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (2,   &r), -3.271093302836352716, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (50,  &r), -37.76583438165180116, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (100, &r), -60.25336482580837088, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (110, &r), -64.2355167606561537,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (111, &r), -64.6268994819519378,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (200, &r), -95.88699147356682665, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_zero_Ai_deriv_e, (2,    &r), -3.248197582179836561, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_deriv_e, (50,   &r), -37.76565910053887108, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_deriv_e, (100,  &r), -60.25329596442479317, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_deriv_e, (110,  &r), -64.23545617243546956, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_deriv_e, (1000, &r), -280.9378080358935071, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (2,    &r), -4.073155089071828216, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (50,   &r), -38.02083574095788210, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (100,  &r), -60.45548887257140819, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (110,  &r), -64.43129648944845060, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (111,  &r), -64.82208737584206093, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (200,  &r), -96.04731050310324450, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_deriv_e, (1000, &r), -281.0315164471118527, TEST_TOL0, GSL_SUCCESS);

  /* Regression vectors for Savannah bug #66808.  The reference values were
   * computed independently from the Maclaurin series in 200-digit
   * arithmetic.  For x < 0 the achievable accuracy is limited by the
   * phase (absolute error ~ eps*|x|^{3/2}), so those use TEST_TOL2. */

  TEST_SF(s, gsl_sf_airy_Ai_e, (-10.0, m, &r),  0.040241238486443190689, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-9.0, m, &r),  -0.022133721547341403674, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-8.0, m, &r),  -0.052705050356386202622, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-7.0, m, &r),   0.18428083525050563728, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-6.0, m, &r),  -0.32914517362982310523, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-5.0, m, &r),   0.35076100902411431979, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-4.0, m, &r),  -0.070265532949289515099, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-3.0, m, &r),  -0.37881429367765807435, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-2.0, m, &r),   0.22740742820168557599, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (-1.0, m, &r),   0.53556088329235211880, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (0.0, m, &r),  0.35502805388781723926, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (1.0, m, &r),  0.13529241631288141552, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (2.0, m, &r),  0.034924130423274379135, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (3.0, m, &r),  0.0065911393574607191395, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (4.0, m, &r),  0.00095156385120480187362, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (5.0, m, &r),  0.00010834442813607441735, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (6.0, m, &r),  9.9476943602528895702e-06, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (7.0, m, &r),  7.4921288639971670808e-07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (8.0, m, &r),  4.6922076160992316256e-08, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (9.0, m, &r),  2.4711684308724898433e-09, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_e, (10.0, m, &r), 1.1047532552898685934e-10, TEST_TOL2, GSL_SUCCESS);

  /* Derivative at the library's own zeros (Savannah bug #66808) */

  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (1, &r), -2.3381074104597670385, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(1), m, &r), 0.70121082272069136249, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (2, &r), -4.087949444130970617, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(2), m, &r), -0.8031113696548639636, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (3, &r), -5.520559828095551059, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(3), m, &r), 0.86520402589415193084, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (4, &r), -6.786708090071758999, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(4), m, &r), -0.9108507370496018030, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (5, &r), -7.944133587120853123, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(5), m, &r), 0.94733570944156776559, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (6, &r), -9.022650853340980380, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(6), m, &r), -0.9779228085694986109, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (7, &r), -10.04017434155808593, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(7), m, &r), 1.0043701226603119685, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (8, &r), -11.00852430373326289, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(8), m, &r), -1.0277386888207861767, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (9, &r), -11.93601556323626252, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(9), m, &r), 1.0487206485881895480, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (10, &r), -12.82877675286575720, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(10), m, &r), -1.0677938591574278346, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (11, &r), -13.69148903521071793, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(11), m, &r), 1.0853028313507000321, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (12, &r), -14.52782995177533498, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(12), m, &r), -1.1015045702774968117, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (13, &r), -15.34075513597799686, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(13), m, &r), 1.1165961779326560845, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (14, &r), -16.13268515694577144, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(14), m, &r), -1.1307323104931878902, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (15, &r), -16.90563399742994263, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(15), m, &r), 1.1440366732735526772, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (16, &r), -17.661300105697057509, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(16), m, &r), -1.1566098491165654638, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (17, &r), -18.401132599207115416, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(17), m, &r), 1.1685347844875248183, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (18, &r), -19.126380474246952144, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(18), m, &r), -1.1798807298701455729, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (19, &r), -19.838129891721499701, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(19), m, &r), 1.1907061311587766528, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Ai_e, (20, &r), -20.537332907677566360, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Ai_deriv_e, (gsl_sf_airy_zero_Ai(20), m, &r), -1.2010607915198232800, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (1, &r), -1.173713222709127925, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(1), m, &r), 0.60195788797623956374, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (2, &r), -3.271093302836352716, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(2), m, &r), -0.7603101414928010888, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (3, &r), -4.830737841662015933, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(3), m, &r), 0.83699101261926109440, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (4, &r), -6.169852128310251260, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(4), m, &r), -0.8894799014265396245, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (5, &r), -7.376762079367763714, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(5), m, &r), 0.92998363856802663447, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (6, &r), -8.491948846509388013, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(6), m, &r), -0.9632344301904237743, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (7, &r), -9.538194379346238887, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(7), m, &r), 0.99158637051766044104, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (8, &r), -10.52991350670535792, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(8), m, &r), -1.0163896592212489424, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (9, &r), -11.47695355127877944, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(9), m, &r), 1.0384942860480093760, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (10, &r), -12.38641713858273875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(10), m, &r), -1.0584718443940233682, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (11, &r), -13.26363952294180555, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(11), m, &r), 1.0767261483165140070, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (12, &r), -14.11275680906865779, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(12), m, &r), -1.0935536233074252884, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (13, &r), -14.93705741215416404, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(13), m, &r), 1.1091786365257356163, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (14, &r), -15.739210351190482771, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(14), m, &r), -1.1237753272862713092, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (15, &r), -16.521419550634379054, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(15), m, &r), 1.1374817081195963948, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (16, &r), -17.285531624581242533, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(16), m, &r), -1.1504091143290550398, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (17, &r), -18.033113287225001572, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(17), m, &r), 1.1626487357883820250, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (18, &r), -18.765508284480081041, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(18), m, &r), -1.1742762531180584491, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (19, &r), -19.483880132989234014, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(19), m, &r), 1.1853552046489041201, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_zero_Bi_e, (20, &r), -20.189244785396202420, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_airy_Bi_deriv_e, (gsl_sf_airy_zero_Bi(20), m, &r), -1.1959394810797524614, TEST_TOL0, GSL_SUCCESS);

  /* #52359: the error estimate used to divide by a series correction that
   * can vanish, which produced err = inf.  Check a sweep of x < -1. */
  {
    double xl;
    for (xl = -1.01; xl > -4.0; xl -= 0.013)
      {
        int status = gsl_sf_airy_Ai_e(xl, m, &r);
        gsl_test_int(status, GSL_SUCCESS, "airy_Ai_e status at %g", xl);
        gsl_test(gsl_isinf(r.err) || gsl_isnan(r.err),
                 "airy_Ai_e finite error at %g", xl);
      }
    {
      int status = gsl_sf_airy_Ai_e(-1.842761151977744, m, &r);
      gsl_test_int(status, GSL_SUCCESS, "airy_Ai_e status at -1.842761151977744");
      gsl_test(gsl_isinf(r.err) || gsl_isnan(r.err),
               "airy_Ai_e finite error at -1.842761151977744");
    }
  }

  /* #52570, #51000: for very large |x| the phase cannot be reduced and the
   * value is indeterminate; the result must stay finite and inside the
   * envelope, and the status must report the loss. */
  {
    const double xhuge[] = { -1.14e34, -1.0e100, -5.643803094122288e102 };
    size_t i;
    for (i = 0; i < sizeof(xhuge) / sizeof(xhuge[0]); ++i)
      {
        double env = 1.0 / sqrt(sqrt(-xhuge[i]));
        int status = gsl_sf_airy_Ai_e(xhuge[i], m, &r);
        gsl_test_int(status, GSL_ELOSS, "airy_Ai_e loss status at %g", xhuge[i]);
        gsl_test(!gsl_finite(r.val), "airy_Ai_e finite value at %g", xhuge[i]);
        gsl_test(fabs(r.val) > env, "airy_Ai_e within envelope at %g", xhuge[i]);

        status = gsl_sf_airy_Bi_e(xhuge[i], m, &r);
        gsl_test_int(status, GSL_ELOSS, "airy_Bi_e loss status at %g", xhuge[i]);
        gsl_test(!gsl_finite(r.val), "airy_Bi_e finite value at %g", xhuge[i]);
        gsl_test(fabs(r.val) > env, "airy_Bi_e within envelope at %g", xhuge[i]);
      }

    {
      double env = sqrt(sqrt(5.643803094122288e102));
      int status = gsl_sf_airy_Ai_deriv_e(-5.643803094122288e102, m, &r);
      gsl_test_int(status, GSL_ELOSS, "airy_Ai_deriv_e loss status");
      gsl_test(!gsl_finite(r.val), "airy_Ai_deriv_e finite value");
      gsl_test(fabs(r.val) > env, "airy_Ai_deriv_e within envelope");
    }

    /* no argument may yield a non-finite result */
    {
      int k;
      for (k = 1; k <= 300; k += 7)
        {
          double x = -pow(10.0, (double) k);
          int status = gsl_sf_airy_Ai_e(x, m, &r);
          gsl_test(!gsl_finite(r.val) || !gsl_finite(r.err),
                   "airy_Ai_e finite at -1e%d", k);
          gsl_test(status != GSL_SUCCESS && status != GSL_ELOSS,
                   "airy_Ai_e status at -1e%d", k);
        }
    }
  }

  return s;
}
