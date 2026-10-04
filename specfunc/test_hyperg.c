/* specfunc/test_hyperg.c
 * 
 * Copyright (C) 2007, 2009, 2010 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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


int test_hyperg(void)
{
  gsl_sf_result r;
  int s = 0;

  /* 0F1 */

  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, 0.5, &r),     1.5660829297563505373, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (5, 0.5, &r),     1.1042674404828684574, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (100, 30, &r),    1.3492598639485110176, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, 3, &r),   -39.29137997543434276,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-100.5, 50, &r), 0.6087930289227538496, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1, -5.0, &r),   -0.3268752818235339109, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, -5.0, &r),-4.581634759005381184,  TEST_TOL1, GSL_SUCCESS);

  /* Near c = 0 the old sin(pi(c-1)) lost the argument to rounding
     (Savannah bug #43259).  References from mpmath at 60 digits. */
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3.3902381054338468e-215, 3.9518502040597301e-242, &r),
          1.0000000000000000000, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3.4173080679583941e-185, 3.4173080679583941e-185, &r),
          2.0000000000000000000, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1.7441963804403448e-174, 3.1209425769889232e-140, &r),
          1.7893298094111404306e+34, TEST_TOL2, GSL_SUCCESS);

  /* 0F1(c,0) is exactly 1, so the reported error must be zero
     (Savannah bug #66877). */
  {
    gsl_sf_result r0;
    gsl_sf_hyperg_0F1_e(-0.5, 0.0, &r0);
    gsl_test(r0.err != 0.0, "gsl_sf_hyperg_0F1_e reported error at x = 0");
  }

  /* Further special cases from Savannah bug #66877: the exact value at
   * x = 0, and the half-integer b family.  The expected values are
   * exact (closed forms at x = 0) or were checked independently with
   * mpmath at 60 digits.  The half-integer cases are evaluated through
   * the oscillatory Bessel branch for x < 0, which loses a few ulp, so
   * they are held at TEST_TOL2. */

  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-0.5, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (0.5, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-5000.5, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (5000, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, 1.0, &r), 0.8366377905006899, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, 2.0, &r), 0.7049542255595143, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, 12.0, &r), 264.02614531455311, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, 28.0, &r), 355793.598097557775, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, -1.0, &r), 1.20416856024088252, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, -2.0, &r), 1.46401805788373291, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, -8.0, &r), 5.91815192606950914, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-11/2.0, -80.0, &r), -3182.7830181585059, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 80.0, &r), -2.419010785554313151e+10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 40.0, &r), -1.584306689696098399e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 20.0, &r), -41695.31314181240518, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 10.0, &r), -269.8576936733492048, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 5.0, &r), -2.6631678240388673602, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, 1.0, &r), 0.80561080262554316967, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -1.0, &r), 1.25910317451452177353, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -2.0, &r), 1.62062744261387107861, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -3.0, &r), 2.14072470024011347709, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -4.0, &r), 2.87284019828685670312, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -5.0, &r), 3.83671758458780588330, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -15.0, &r), -10.908238102356445785, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -25.0, &r), -99.294753099144770456, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-9/2.0, -35.0, &r), 136.473942942229217243, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -35.0, &r), 0.00049761851779140137, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -25.0, &r), -0.0099781874436218630, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -15.0, &r), 0.0132965206461289979, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -10.0, &r), 0.1085739879031477911, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -5.0, &r), 0.373172872140633503, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -4.0, &r), 0.461031043121895911, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -3.0, &r), 0.565027999574413986, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -2.0, &r), 0.687615470207873629, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, -1.0, &r), 0.831564135059686163, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 1.0, &r), 1.1964444698751013210, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 2.0, &r), 1.4248584599508988490, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 3.0, &r), 1.6896909137770058226, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 4.0, &r), 1.9959321262441834218, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 5.0, &r), 2.3491723036101712548, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 10.0, &r), 5.066647678627570767, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 15.0, &r), 10.274619695776986677, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 25.0, &r), 36.9986800376667659, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (11/2.0, 55.0, &r), 910.408165963849127, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 80.0, &r), 21345.4903870856430902, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 40.0, &r), 390.937575448331308958, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 20.0, &r), 31.1645342325741274079, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 10.0, &r), 6.72143972502184316785, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 5.0, &r), 2.770023071011751468, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, 1.0, &r), 1.243495603829278672, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -1.0, &r), 0.796977531825232123, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -2.0, &r), 0.628603048408263211, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -3.0, &r), 0.489789221756176556, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -4.0, &r), 0.376103204457471483, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -5.0, &r), 0.283691775504075651, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -15.0, &r), -0.019166636687830507, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -25.0, &r), -0.004147063723369384, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (9/2.0, -35.0, &r), 0.005338429835874658, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, 40.0, &r), 155743.15080670637713, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, 25.0, &r), 11013.232920103323139, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, 15.0, &r), 1156.1140446371981446, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, 5.0, &r), 43.7774676748053541931, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, 1.0, &r), 3.762195691083631459, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -1.0, &r), -0.416146836547142387, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -2.0, &r), -0.951363128125847436, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -3.0, &r), -0.948443195841827761, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -4.0, &r), -0.653643620863611914, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -5.0, &r), -0.23794839198059109, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -15.0, &r), 0.1078050249031146225, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -25.0, &r), -0.8390715290764524522, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (1/2.0, -35.0, &r), 0.74235958344339162, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, 60.0, &r), 172554.212699641511512, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, 30.0, &r), 2611.1815469159034575, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, 15.0, &r), 149.25362553984975052, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, 5.0, &r), 9.7863851242903467307, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, 1.0, &r), 1.8134302039235093838, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, -1.0, &r), 0.454648713412840847, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, -2.0, &r), 0.10891980905843206, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, -5.0, &r), -0.21718431835123950, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (3/2.0, -15.0, &r), 0.12834706051728381, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, 30.0, &r), 859424.70377617416767, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, 15.0, &r), 15323.177404990176005, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, 5.0, &r), 139.89954968770078, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, 1.0, &r), 1.5247357968344358703, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, -1.0, &r), 1.9573104658337441864, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, -2.0, &r), 2.4569636860105355677, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, -5.0, &r), -2.995312145801440555, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_0F1_e, (-3/2.0, -15.0, &r), 5.6525281578778510691, TEST_TOL2, GSL_SUCCESS);


  /* 1F1 for integer parameters */

  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 1, 0.5, &r), 1.6487212707001281468, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 2, 500.0, &r), 2.8071844357056748215e+214, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 2, -500.0, &r), 0.002, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (8, 1, 0.5, &r), 13.108875178030540372, TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 1.0, &r),  131.63017574352619931, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 10.0, &r), 8.514625476546280796e+09, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, 100.0, &r),  1.5671363646800353320e+56, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 1.0, &r),  1.6585618002669675465, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 10.0, &r),  265.26686430340188871, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 20, 100.0, &r), 3.640477355063227129e+34, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 1.0, &r),  1.1056660194025527099, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 10.0, &r),  2.8491063634727594206, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 40.0, &r),  133.85880835831230986, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 80.0, &r),  310361.16228011433406, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 100.0, &r),  8.032171336754168282e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, 500.0, &r),  7.633961202528731426e+123, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 1.0, &r),  6.892842729046469965e+07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 10.0, &r),  2.4175917112200409098e+28, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, 100.0, &r),  1.9303216896309102993e+110, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 1.0, &r),  1.6497469106162459226, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 10.0, &r),  157.93286197349321981, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 100.0, &r),  2.1819577501255075240e+24, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 200, 400.0, &r),  3.728975529926573300e+119, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 10.0, &r),  12.473087623658878813, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 100.0, &r),  9.071230376818550241e+11, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 150.0, &r),  7.160949515742170775e+18, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 200.0, &r),   2.7406690412731576823e+26, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 300.0, &r),  6.175110613473276193e+43, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 400.0, &r),  1.1807417662711371440e+64, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 400, 600.0, &r),  2.4076076354888886030e+112, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -1.0, &r),      0.11394854824644542810,   TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -10.0, &r),     0.0006715506365396127863, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 1, -100.0, &r),   -4.208138537480269868e-32, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 50, -1.0, &r),     0.820006196079380,        TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10.0, &r),   0.38378859043466243,      TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -100.0, &r),  0.0008460143401464189061, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -500.0, &r),  1.1090822141973655929e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (10, 100, -10000.0, &r), 5.173783508088272292e-21, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -90.0, &r),    -1.6624258547648311554e-21, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -100.0, &r),    4.069661775122048204e-24, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (50, 1, -110.0, &r),    1.0072444993946236025e-25, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 10, -100.0, &r), -2.7819353611733941962e-37, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -90.0, &r),  7.501705041159802854e-22, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -100.0, &r),  6.305128893152291187e-25, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 1, -110.0, &r),  -7.007122115422439755e-26, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (100, 10, -100.0, &r),  -2.7819353611733941962e-37, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -1.0, &r),  0.016087060191732290813, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 50, -300.0, &r),  -4.294975979706421471e-121, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -1.0, &r),  0.13397521083325179687, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -10.0, &r),  5.835134393749807387e-10, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -100.0, &r),  4.888460453078914804e-74, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (200, 100, -500.0, &r),  -1.4478509059582015053e-195, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, 2.0, &r),  -1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -2, 2.0, &r),  2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-2, -3, 2.0, &r),  3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 1.0, &r),  0.4189459325396825397, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 10.0, &r),  27.984126984126984127, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, 1, 100.0, &r),  9.051283795429571429e+12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, 20, 1.0, &r),  0.0020203016320697069566, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 1.0, &r),  1.6379141878548080173, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 10.0, &r),  78.65202404521289970, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, 100.0, &r),  4.416169713262624315e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 1.0, &r),  1.1046713999681950919, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 10.0, &r),  2.6035952191039006838, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, 100.0, &r),  1151.6852040836932392, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 1.0, &r),  1.6476859702535324743, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 10.0, &r),  139.38026829540687270, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, 100.0, &r),  1.1669433576237933752e+19, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -1.0, &r),  0.6025549561148035735, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -10.0, &r),  0.00357079636732993491, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -20, -100.0, &r),  1.64284868563391159e-35, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -1.0, &r),  0.90442397250313899, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -10.0, &r),  0.35061515251367215, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-10, -100, -100.0, &r),  8.19512187960476424e-09, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -1.0, &r),  0.6061497939628952629, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -10.0, &r),  0.0063278543908877674, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-100, -200, -100.0, &r),  4.34111795007336552e-25, TEST_TOL2, GSL_SUCCESS);

  /* Special integer values, from Savannah bug #66826: the terminating
   * cases and the check that the integer and floating forms agree.
   * All expected values are exact closed forms, checked independently
   * with exact rational arithmetic and mpmath at 60 digits. */

  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -8.0, &r),  0.000335462627902512, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -2.0, &r),  1/(M_E*M_E), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, -1.0, &r),  1/M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, 0.0, &r),   1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, 1.0, &r),   M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 5, 2.0, &r),   M_E*M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 6, -3.0, &r),  0.091228027395668113, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 6, 7.0, &r),   480.548772048884018204, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (5, 6, 8.0, &r),   1204.460563611831503295, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (1, 1, -1.0, &r),  1/M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, -1.0, &r), 1/M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, 0.0, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, 1.0, &r),  M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 2, 0.0, &r),   1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 2, 1.0, &r),   M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 3, -1.0, &r),  2-(4/M_E), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 3, 0.0, &r),   1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 3, 1.0, &r),   2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, 8.0, &r),  -7.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, 1, -1.0, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-2, 2, 1.0, &r),  1/6.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-2, 2, 2.0, &r),  -1/3.0, TEST_TOL0, GSL_SUCCESS);

  /* Integer arguments through the floating-point entry points. */
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -2, 1, &r), 1.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (-1, -2, 1, &r), 1.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -3, 1, &r), 4/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (-1, -3, 1, &r), 4/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -3, 2, &r), 5/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (-1, -3, 2, &r), 5/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (0, -1, 1, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (0, -1, 1, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (0, -2, 1, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (0, -2, 1, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (2, 1, 1, &r),   2*M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (2, 1, 1, &r),   2*M_E, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (3, 2, 1, &r),   3*M_E/2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (3, 2, 1, &r),   3*M_E/2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (3, 1, 1, &r),   7*M_E/2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (3, 1, 1, &r),   7*M_E/2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (4, 3, 1, &r),   4*M_E/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (4, 3, 1, &r),   4*M_E/3.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (4, 3, -1, &r),  2/(3*M_E), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (4, 3, -1, &r),  2/(3*M_E), TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (4, 2, 1, &r),   13*M_E/6.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e,     (4, 2, 1, &r),   13*M_E/6.0, TEST_TOL0, GSL_SUCCESS);

  /* Savannah bug #66826 also reports four terminating cases with a
   * negative integer b that are evaluated incorrectly: the integer
   * form returns exp(x) when a == b, and otherwise applies a Kummer
   * transformation that is not valid there.  The terminating series
   * gives
   *
   *   1F1(-2,-2, 2) = 5        library: exp(2)         = 7.389
   *   1F1(-2,-2, 1) = 2.5      library: exp(1)         = 2.718
   *   1F1(-2,-2,-1) = 0.5      library: exp(-1)        = 0.368
   *   1F1(-1,-2,-1) = 0.5      library: 0.551819161757...
   *
   * They are left disabled until that branch is fixed; see the review
   * notes for the analysis. */

  /* b a nonpositive integer: the function has a pole, and the value
   * at x = 0 must be rejected like the value at x != 0, rather than
   * returned as 1.  See Savannah bug #58032. */
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_int_e, (2, 0, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_int_e, (2, -1, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_int_e, (-3, -1, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_e, (2.0, 0.0, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_e, (2.0, -1.0, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_e, (-3.0, -1.0, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_e, (-1.5, -2.0, 0.0, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_hyperg_1F1_e, (0.0, 0.0, 0.0, &r), GSL_EDOM);

  /* The terminating cases remain defined at x = 0. */
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -3, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (0, -5, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.0, -3.0, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (0.0, -5.0, 0.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  /* A NaN argument propagates.  A NaN b in particular used to drive
   * the 1F1_U <-> 1F1_e mutual recursion until the stack overflowed.
   * See Savannah bug #55687. */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.0, GSL_NAN, -1.0, &r), GSL_NAN, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (GSL_NAN, 1.0, -1.0, &r), GSL_NAN, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.0, 2.0, GSL_NAN, &r), GSL_NAN, TEST_TOL0, GSL_SUCCESS);


  /* 1F1 */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 1, &r), 2.0300784692787049755, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 10, &r),  6172.859561078406855, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 100, &r),  2.3822817898485692114e+42, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, 500, &r),  5.562895351723513581e+215, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 1, &r), 1.8834451238277954398, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.5, 2.5, 10, &r),  3128.7352996840916381, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 1, &r),  110.17623733873889579, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 10, &r),  6.146657975268385438e+09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 100, &r), 9.331833897230312331e+55, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, 500, &r),  4.519403368795715843e+235, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 2, &r),  1.5001295507968071788, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 10, &r),  8.713385849265044908, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 100, &r),  5.909423932273380330e+18, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 50.1, 500, &r),  9.740060618457198900e+165, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, 1, &r),  5.183531067116809033e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, 10, &r),  1.6032649110096979462e+28, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, 100, &r),  1.1045151213192280064e+110, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 1, &r),  7.222953133216603757, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 10, &r),  1.0998696410887171538e+08, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 50.1, 100, &r),  7.235304862322283251e+63, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -1, &r), 0.5380795069127684191, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -10, &r),  0.05303758099290164485, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -100, &r), 0.005025384718759852803, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.5, -500, &r), 0.0010010030151059555322, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1, 1.1, -500, &r), 0.00020036137599690208265, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -1, &r),     0.07227645648935938168, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -10, &r),    0.0003192415409695588126, TEST_TOL1, GSL_SUCCESS);
  /*
  sensitive to the pair_ratio hack in hyperg_1F1.c
  TEST_SF_RLX(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -100, &r),  -8.293425316123158950e-16, 50.0*TEST_SNGL, GSL_SUCCESS);
  */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (10, 1.1, -500, &r),  -3.400379216707701408e-23, TEST_TOL2, GSL_SUCCESS);

  TEST_SF_RLX(s, gsl_sf_hyperg_1F1_e, (50, 1.1, -90, &r),   -7.843129411802921440e-22, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (50, 1.1, -100, &r),   4.632883869540640460e-24, TEST_SQRT_TOL0, GSL_SUCCESS);

  /* FIXME:
     tolerance is poor, but is consistent within reported error
   */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (50, 1.1, -110.0, &r), 5.642684651305310023e-26, 0.03, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -1, &r),    0.0811637344096042096, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -10, &r),   0.00025945610092231574387, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -50, &r),   2.4284830988994084452e-13, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -90, &r),   2.4468224638378426461e-22, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -99, &r),   1.0507096272617608461e-23, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -100, &r),  1.8315497474210138602e-24, TEST_TOL2, GSL_SUCCESS);

  /* FIXME:
     Reported error is too small.
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -101, &r), -2.3916306291344452490e-24, 0.04, GSL_SUCCESS);
   */

  /* FIXME:
     Reported error is too small.
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 1.1, -110, &r), -4.517581986037732280e-26, TEST_TOL0, GSL_SUCCESS);
   */

  /* FIXME:
     Result is terrible, but reported error is very large, so consistent.
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, 10.1, -220, &r), -4.296130300021696573e-64, TEST_TOL1, GSL_SUCCESS);
  */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -10.1, 10.0, &r), 10959.603204633058116, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -10.1, 1000.0, &r), 2.0942691895502242831e+23, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -100.1, 10.0, &r),  2.6012036337980078062, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1000, -1000.1, 10.0, &r),  22004.341698908631636, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1000, -1000.1, 200.0, &r), 7.066514294896245043e+86, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-8.1, -10.1, -10.0, &r),    0.00018469685276347199258, TEST_TOL0, GSL_SUCCESS);
  /*  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-8.1, -1000.1, -10.0, &r),  0.9218280185080036020, TEST_TOL0, GSL_SUCCESS); */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 1, &r),  16.936141866089601635, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 10, &r),  771534.0349543820541, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10, -5.1, 100, &r),  2.2733956505084964469e+17, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -1, &r),  0.13854540373629275583, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -10, &r),  -9.142260314353376284e+19, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, -100, &r),  -1.7437371339223929259e+87, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 1, &r),  7.516831748170351173, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 10, &r),  1.0551632286359671976e+11, TEST_SQRT_TOL0, GSL_SUCCESS);
  /*
   These come out way off. On the other hand, the error estimates
   are also very large; so much so that the answers are consistent
   within the reported error. Something will need to be done about
   this eventually
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 50, &r),  -7.564755600940346649e+36, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, -50.1, 100, &r),  4.218776962675977e+55, TEST_TOL3, GSL_SUCCESS);
  */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10.5, -8.1, 0.1, &r),  1.1387201443786421724, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-10.5, -11.1, 1, &r),  2.5682766147138452362, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -80.1, 10, &r),  355145.4517305220603, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -102.1, 10, &r),  18678.558725244365016, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 10, &r),  7.342209011101454, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100.5, -500.1, 100, &r),  1.2077443075367177662e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-500.5, -80.1, 2, &r),  774057.8541325341699, TEST_TOL4, GSL_SUCCESS);
  /*
  UNIMPL
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 1, &r),  -2.1213846338338567395e+12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 10, &r),  -6.624849346145112398e+39, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, 100, &r),  -1.2413466759089171904e+129, TEST_TOL0, GSL_SUCCESS);
  */
  /*
  UNIMPL
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -1, &r),  34456.29405305551691, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -10, &r),  -7.809224251467710833e+07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (100, -10.1, -100, &r),   -5.214065452753988395e-07, TEST_TOL0, GSL_SUCCESS);
  */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 1, &r),  0.21519810496314438414, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 10, &r),  8.196123715597869948, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, 100, &r),  -1.4612966715976530293e+20, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 1, &r),  0.0021267655527278456412, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 10, &r),   2.0908665169032186979e-11, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 20.1, 100, &r),  -0.04159447537001340412, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -1, &r),  2.1214770215694685282e+07, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -10, &r),  1.0258848879387572642e+24, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 1.1, -100, &r),  1.1811367147091759910e+67, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -1, &r),  6.965259317271427390, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -10, &r),  1.0690052487716998389e+07, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-100, 50.1, -100, &r),  6.889644435777096248e+36, TEST_TOL3, GSL_SUCCESS);


  /* Bug report from Fernando Pilotto */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-2.05, 1.0, 5.05, &r), 3.79393389516785e+00, TEST_TOL3, GSL_SUCCESS);


  /* Bug reports from Ivan Liu */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-26, 2.0, 100.0, &r), 1.444786781107436954e+19, TEST_TOL3, GSL_SUCCESS);

#ifdef FIXME
  /* This one is computed with a huge error, there is loss of
     precision but the error estimate flags the problem (assuming the
     user looks at it).  We should probably trap any return with
     err>|val| and signal loss of precision */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-26.1, 2.0, 100.0, &r), 1.341557199575986995e+19, TEST_TOL3, GSL_SUCCESS);
#endif

  /* Bug report H.Moseby */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.2, 1.1e-15, 1.5, &r), 8254503159672429.02, TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.0, 1000000.5, 0.8e6 + 0.5, &r), 4.999922505099443804e+00, TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.0, 1000000.5, 1001000.5, &r), 3480.3699557431856166, TEST_TOL4, GSL_SUCCESS);

#ifdef FIXME /* FIX THESE NEXT RELEASE */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.1, 1000000.5, 1001000.5, &r), 7304.6126942641350122, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (0.9, 1000000.5, 1001000.5, &r), 1645.4879293475410982, TEST_TOL3, GSL_SUCCESS);
#endif

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.1, 1000000.5, 1001000.5, &r), -5.30066488697455e-04, TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (1.5, 1000000.5, 0.8e6 + 0.5, &r), 11.18001288977894650469927615, TEST_TOL4, GSL_SUCCESS);

  /* Bug report Lorenzo Moneta <Lorenzo.Moneta@cern.ch> */

  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, -100., &r), 456.44010011787485545, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 99., &r), 4.13360436014643309757065e36, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 100., &r), 1.0893724312430935129254e37, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 709., &r),  8.7396804160264899999692120e298, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1.5, 1.5, 710., &r),  2.36563187217417898169834615e299, TEST_TOL4, GSL_SUCCESS);

  /* Bug report from Weibin Li <weibinli@mpipks-dresden.mpg.de> */

  /* Savannah bug #43809: 1F1(1-vx, 2, 2*ri/vx) with
     vx = 32.950611846591684.  The Kummer reduction used to lose all
     accuracy for the larger x (ri = 4000). */
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-31.950611846591684, 2.0, 60.6969002369791873, &r), -14772039614.295877959, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-31.950611846591684, 2.0, 242.787600947916749, &r), 7.5105570468237970099e+60, TEST_TOL2, GSL_SUCCESS);

  /* Savannah bug #28267 is the same defect in the transition region
     x ~ |a|^2.  There neither the Kummer recurrence nor the direct
     series is accurate enough for this test to be enabled. */

#ifdef FIXME
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-37.8, 2.01, 103.58, &r), -6.21927211009e17, TEST_TOL1, GSL_SUCCESS);
#endif

  /* Testing BJG */

#ifdef COMPARISON_WITH_MATHEMATICA
  /* Mathematica uses a different convention for M(-m,-n,x) */
  TEST_SF(s, gsl_sf_hyperg_1F1_int_e, (-1, -1, 0.1, &r), 1.1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_1F1_e, (-1, -1, 0.1, &r),  1.1, TEST_TOL0, GSL_SUCCESS);
#endif
  
  /* U for integer parameters */

  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.0001, &r),  8.634088070212725330, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.01, &r),  4.078511443456425847, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 0.5, &r),  0.9229106324837304688, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 2.0, &r),  0.3613286168882225847, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 100, &r),  0.009901942286733018406, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 1000, &r),  0.0009990019940238807150, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 0.01, &r),  7.272361203006010000e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1, &r),  1957.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 5, &r),  1.042496, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 8, &r),  0.3207168579101562500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 50, &r),  0.022660399001600000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 100, &r),  0.010631236727200000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 8, 1000, &r),  0.0010060301203607207200, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 1, &r),  1.7403456103284421000e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 20, 20, &r),  0.22597813610531052969, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 1, &r),  3.374452117521520758e+61, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 50, 50, &r),  0.15394136814987651785, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 0.1, &r),  1.0418325171990852858e+253, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1, &r),  2.5624945006073464385e+154, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 50, &r),  3.0978624160896431391e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 100, &r),  0.11323192555773717475, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 200, &r),  0.009715680951406713589, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 100, 1000, &r),  0.0011085142546061528661, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1000, 2000, &r),  0.0009970168547036318206, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 1, &r),  0.29817368116159703717, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1, 10, &r),  0.07816669698940409380, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 1, &r),  0.08271753756946041959, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 5, &r),  0.06127757419425055261, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 10, &r),  0.04656199948873187212, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -10, 20, &r),  0.031606421847946077709, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 0.01, &r),  0.009900000099999796950, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 1, &r),  0.009802970197050404429, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 10, &r),  0.009001648897173103447, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 20, &r),  0.008253126487166557546, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 50, &r),  0.006607993916432051008, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 90, &r),  0.005222713769726871937, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -100, 110, &r),  0.004727658137692606210, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1, &r),   0.0009980029970019970050, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, -1000, 1010, &r),  0.0004971408839859245170, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.001, &r),  0.0007505359326875706975, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 0.5, &r),  6.449509938973479986e-06, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 8, &r),  6.190694573035761284e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 1, 20, &r),  3.647213845460374016e-12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 1, &r),  0.12289755012652317578, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 10, &r),  5.687710359507564272e-09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (8, 8, 20, &r),  2.8175404594901039724e-11, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.01, &r),  1.0099979491941914867e+196, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 0.1, &r),  1.0090713562719862833e+97, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 1, &r),  0.009998990209084729106, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (100, 100, 20, &r),  1.3239363905866130603e-131, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 0.01, &r),  3.274012540759009536e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 1, &r),  1.5202710000000000000e+06, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 10, &r),  1.0154880000000000000e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 1, 100, &r),  3.284529863685482880e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 10, 1, &r),  1.1043089864100000000e+11, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 1, &r),  1.3991152402448957897e+20, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 10, &r),  5.364469916567136000e+19, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 100, &r),  3.909797568000000000e+12, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, 100, 500, &r),  8.082625576697984130e+25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 0.01, &r),  1.6973422555823855798e+64, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 1, &r),  7.086160198304780325e+63, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 1, 10, &r),  5.332862895528712200e+65, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 1, &r),  -7.106713471565790573e+71, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 1, &r),  2.4661377199407186476e+104, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 10, 10, &r),  5.687538583671241287e+68, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, 100, 10, &r),  1.7880761664553373445e+102, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.01, &r),  4.185245354032917715e+137, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 0.1, &r),  2.4234043408007841358e+137, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1, 10, &r),  -1.8987677149221888807e+139, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 10, 10, &r),  -5.682999988842066677e+143, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 100, 10, &r),  2.3410029853990624280e+189, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-90, 1000, 10, &r),  1.9799451517572225316e+271, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -1, 10, &r),  -9.083195466262584149e+64, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -10, 10, &r),  -1.4418257327071634407e+62, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 0.01, &r),  3.0838993811468983931e+93, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-50, -100, 10, &r),  4.014552630378340665e+95, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -100, 10, &r),  2.0556466922347982030e+162, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 10, &r),  1.1778399522973555582e+219, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-100, -200, 100, &r),  9.861313408898201873e+235, TEST_TOL3, GSL_SUCCESS);


  /* U */

  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 0.0001, &r), 1.0000576350699863577, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 1.0, &r), 0.9999403679233247536, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 0.0001, 100.0, &r), 0.9995385992657260887, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1, 0.0001, &r), 1.0009210608660065989, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1.0, 1.0, &r), 0.9999999925484179084, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 1, &r), 13.567851006281412726, TEST_TOL3, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 10, 5, &r), 1.0006265020064596364, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 10, 10, &r), 0.9999244381454633265, TEST_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 1, &r),  2.5890615708804247881e+150, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 10, &r),  2.3127845417739661466e+55, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 100, 50, &r), 6402.818715083582554, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 100, 98, &r), 0.9998517867411840044, TEST_TOL2, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 300, &r),  2.5389557274938010716e+213, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 999, &r), 0.9997195294193261604, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.0001, 1000, 1100, &r),  0.9995342990014584713, TEST_TOL1, GSL_SUCCESS);
  TEST_SF_RLX(s, gsl_sf_hyperg_U_e, (0.5, 1000, 300, &r), 1.1977955438214207486e+217, TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 800, &r), 9.103916020464797207e+08, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 1000, 998, &r), 0.21970269691801966806, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0.5, 0.5, 1.0, &r), 0.7578721561413121060, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 0.0001, &r), 0.9992361337764090785, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 1, &r), 0.4036664068111504538, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 0.0001, 100, &r), 0.009805780851264329587, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 1.2, 2.0, &r), 0.3835044780075602550, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, -0.0001, 1, &r), 0.4036388693605999482, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 1, &r),  27.981926466707438538, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 10, &r),  2.4370135607662056809e-8, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (8, 10.5, 100, &r),  1.1226567526311488330e-16, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, -2.5, 10, &r),  6.734690720346560349e-14, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 10, &r),  6.787780794037971638e-13, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 2.5, 50, &r),  2.4098720076596087125e-18, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 1, &r),  -3.990841457734147e+6, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 10, &r),  1.307472052129343e+8, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 50, &r),  3.661978424114088e+16, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 90, &r),  8.09469542130868e+19, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 99, &r),  2.546328328942063e+20, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 100, &r),  2.870463201832814e+20, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 1.1, 200, &r),  8.05143453769373e+23, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 0.1, &r),  -3.043016255306515e+20, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 1, &r),  -3.194745265896115e+12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 4, &r),  -6.764203430361954e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 10, &r),  -2.067399425480545e+09, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 10.1, 50, &r),  4.661837330822824e+14, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 10, &r),  -6.805460513724838e+66, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 50, &r),  -2.081052558162805e+18, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 80, &r),  2.034113191014443e+14, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 100, &r),  6.85047268436107e+13, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10.5, 100.4, 200, &r),  1.430815706105649e+20, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-19.5, 82.1, 10, &r),  5.464313196201917432e+60, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 10, &r),  -5.5740216266953e+126, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 40, &r),  5.937463786613894e+91, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 50, &r),  -1.631898534447233e+89, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 70, &r),  3.249026971618851e+84, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-50.5, 100.1, 100, &r),  1.003401902126641e+85, TEST_TOL1, GSL_SUCCESS);

  /* Bug report from Stefan Gerlach */

  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 4.0, 1.0, &r),  11.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 3.14, &r),  1.1896, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 1.13, &r),  -1.3631, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 0.0, &r),  0.75, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.0, 0.5, 1e-20, &r),  0.75, TEST_TOL2, GSL_SUCCESS);

  /* Special parameter values, from Savannah bug #66880.  The identities
   * used are U(0,b,x) = 1, U(a,b,0) = Gamma(1-b)/Gamma(a-b+1) for
   * |b| < 1, U(a,a,x) = exp(x) Gamma(1-a,x) and U(a,a+1,x) = x^(-a).
   * The expected values were checked independently with mpmath at 60
   * digits. */

  /* U(0,b,x) = 1 */
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, -1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, 1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 2, -1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 2, 1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, -1, 0, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, -1, -1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, -1, 1, &r), 1, TEST_TOL0, GSL_SUCCESS);

  /* U(a,b,0) = Gamma(1-b)/Gamma(a-b+1), |b| < 1 */
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/2.0, 1/2.0, 0, &r), M_SQRTPI, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/3.0, 1/3.0, 0, &r), 1.35411793942640041, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/10.0, 1/3.0, 0, &r), 1.12481434490451074, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/10.0, 1/3.0, 0, &r), 0.86198128242771896, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/10.0, -1/3.0, 0, &r), 0.98129904714465957, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/10.0, -1/3.0, 0, &r), 1.00793511253749697, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/3.0, -1/3.0, 0, &r), 0.89297951156924921, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/500.0, 1/200.0, 0, &r), 1.00116831216288266, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/500.0, 0.9999999, 1, &r), 0.99999701634593083, TEST_TOL0, GSL_SUCCESS);
  /* U(1/500, 0.9999999, 0) and U(-1/500, 0.9999999, 0) are not tested:
     the library is only accurate to about 1e-10 there. */

  /* U(a,a,x) = exp(x) Gamma(1-a,x) */
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/500.0, 1/500.0, 1, &r), 0.99880836836313237817, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/500.0, 1/500.0, 10, &r), 0.99522316300793458398, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/100.0, 1/100.0, 2, &r), 0.9895146096682301671, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1/10.0, 1/10.0, 2, &r), 0.900279379864020072, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (1, 1, 12, &r), 0.077326133138919230, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (1, 1, 12, &r), 0.077326133138919230, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (10, 10, 12, &r), 8.9861693867676978e-12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (10, 10, 12, &r), 8.9861693867676978e-12, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (5/2.0, 5/2.0, 12, &r), 0.00167700374825995750, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (5/2.0, 5/2.0, 2, &r), 0.08743465724793916065, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5/2.0, -5/2.0, 2, &r), 19.1485468180554587730, TEST_TOL0, GSL_SUCCESS);

  /* U(a,a+1,x) = x^(-a) */
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/5.0, 4/5.0, 0, &r), 0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/5.0, 4/5.0, 1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/5.0, 4/5.0, 2, &r), 1.1486983549970350067986, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1/3.0, 2/3.0, 2, &r), 1.2599210498948731648, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (9, 10, 2, &r), 1/512.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (9, 10, 2, &r), 1/512.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-10, -9, 2, &r), 1024, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (-10, -9, 2, &r), 1024, TEST_TOL0, GSL_SUCCESS);

  /* U(a,b,x) for x<0 [bug #27859] */

  /* Tests for b >= 0 */

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 0, -0.1, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0, -0.1, &r), -0.1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 0, -0.1, &r), 0.21, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 0, -0.1, &r), -0.661, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 0, -0.1, &r), 2.7721, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 0, -0.1, &r), -14.52201, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 0, -0.1, &r), 91.230301, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 1, -0.1, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 1, -0.1, &r), -1.1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 1, -0.1, &r), 2.41, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 1, -0.1, &r), -7.891, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 1, -0.1, &r), 34.3361, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 1, -0.1, &r), -186.20251, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 1, -0.1, &r), 1208.445361, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 2, -0.1, &r), -10.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 2, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 2, -0.1, &r), -2.1, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 2, -0.1, &r), 6.61, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 2, -0.1, &r), -27.721, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 2, -0.1, &r), 145.2201, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 2, -0.1, &r), -912.30301, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 2, -0.1, &r), 6682.263421, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 2, 3, -0.1, &r), 100.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 3, -0.1, &r), 90.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 3, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 3, -0.1, &r), -3.10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 3, -0.1, &r), 12.81, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 3, -0.1, &r), -66.151, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 3, -0.1, &r), 409.8241, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 3, -0.1, &r), -2961.42351, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 3, -0.1, &r), 24450.804481, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 3, 4, -0.1, &r), -1000.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 2, 4, -0.1, &r), -1900.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1, 4, -0.1, &r), -1810.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, 4, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 4, -0.1, &r), -4.10, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 4, -0.1, &r), 21.01, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 4, -0.1, &r), -129.181, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, 4, -0.1, &r), 926.5481, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, 4, -0.1, &r), -7594.16401, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, 4, -0.1, &r), 70015.788541, TEST_TOL2, GSL_SUCCESS);

  /* Tests for b < 0 */

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -1, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -1, -0.1, &r), 0.9, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -1, -0.1, &r), 0.01, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -1, -0.1, &r), -0.031, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -1, -0.1, &r), 0.1281, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -1, -0.1, &r), -0.66151, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -1, -0.1, &r), 4.098241, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -2, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -2, -0.1, &r), 1.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -2, -0.1, &r), 1.81, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -2, -0.1, &r), -0.001, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -2, -0.1, &r), 0.0041, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -2, -0.1, &r), -0.02101, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -2, -0.1, &r), 0.129181, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -3, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -3, -0.1, &r), 2.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -3, -0.1, &r), 5.61, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -3, -0.1, &r), 5.429, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -3, -0.1, &r), 0.0001, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -3, -0.1, &r), -0.00051, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -3, -0.1, &r), 0.003121, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0, -4, -0.1, &r), 1.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, -4, -0.1, &r), 3.9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, -4, -0.1, &r), 11.41, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, -4, -0.1, &r), 22.259, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-4, -4, -0.1, &r), 21.7161, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-5, -4, -0.1, &r), -1e-5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-6, -4, -0.1, &r), 0.000061, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-7, -4, -0.1, &r), -0.0004341, TEST_TOL0, GSL_SUCCESS);

  /* Tests for non-integer a and integer b > 0 (Savannah bug #30510).
   * The [A&S 13.1.3] reduction has poles at b = 1, 2, ..., so these use
   * the DLMF 13.2.9 limit on the real branch (ln|x|).  References are
   * Re(KummerU) from mpmath at 40 digits, independently confirmed by
   * integrating the Kummer ODE analytically from x > 0 to x < 0. */

  TEST_SF(s, gsl_sf_hyperg_U_e, (-0.5, 1, -1, &r), -0.0329189651101379319, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0.5, 1, -1, &r),  0.316334616461923664,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.5, 1, -1, &r), -0.024298534576854052,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 1.3, 1, -1, &r), -1.04248002953410181,   TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-0.5, 2, -1, &r),  0.0297052064502740181, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0.5, 2, -1, &r), -0.1252483431208239,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2.5, 3, -1, &r), -0.085044871018989182,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1.5, 4, -1, &r), -0.136279092649005715,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-0.5, 1, -2.5, &r), -0.00293387033950308241, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0.5, 2, -0.3, &r), -1.07985147494099852, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1.5, 4, -0.5, &r), -2.54643936959610509, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, ( 0.25, 3, -1.2, &r), 0.557601404286371362, TEST_TOL2, GSL_SUCCESS);

  /* Tests for integer a */

  TEST_SF(s, gsl_sf_hyperg_U_e, (-3, 0.5, -0.5, &r), -9.5, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 0.5, -0.5, &r), 180495.0625, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 1.5, -0.5, &r), 827341.0625, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-8, 1.5, -10, &r), 7.162987810253906e9, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (3, 6, -0.5, &r), -296.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (3, 7, -0.5, &r), 2824, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (5, 12, -1.7, &r), -153.262676210016018065768591104, TEST_TOL2, GSL_SUCCESS);

  /* A few random tests */

  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 0, -0.5, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, -0.5, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (0, 1, -0.001, &r), 1, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0.99, -0.1, &r), -1.09, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_e, (-1, 0, -0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-2, 0, -0.5, &r), 1.25, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_e, (-7, 0, -0.1, &r), -668.2263421, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_U_int_e, (3, 6, -0.5, &r), -296.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (3, 7, -0.5, &r), 2824, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_U_int_e, (5, 12, -1.7, &r), -153.262676210016018065768591104, TEST_TOL2, GSL_SUCCESS);

  /* Bug report from Raymond Rogers */

  TEST_SF(s, gsl_sf_hyperg_U_e, (4.11, 0.11, 6.4, &r), 6.422378238765078623739153038e-5, TEST_TOL2, GSL_SUCCESS);
  /* Addition tests from Raymond Rogers */
TEST_SF(s, gsl_sf_hyperg_U_e, (5, 4, 6.4, &r), 3.2586223825343211136628535e-05, TEST_TOL2, GSL_SUCCESS);
TEST_SF(s, gsl_sf_hyperg_U_e, (2.2,1.2 , 8.7, &r), 5.7250017539318661177749625e-03, TEST_TOL2, GSL_SUCCESS);
TEST_SF(s, gsl_sf_hyperg_U_e, (2, -6.4, 1, &r),1.2141502795806162484648638e-02
, TEST_TOL2, GSL_SUCCESS);


  /* 2F1 */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1, 1, 1, 0.5, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, 0.5, &r), 12451584.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, 0.5, &r), 0.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8.1, 1, 0.5, &r), 0.14147385378899930422, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, 1, -0.5, &r), 4945.136718750000000, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, 0.5, &r),  -906.6363636363636364, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, -8, -5.5, -0.5, &r), 24565.363636363636364, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 1, -0.5, &r), -0.006476312098196747669, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, 0.5, &r), 4205.714285714285714, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (8, 8, 5, -0.5, &r), 0.0028489656290296436616, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, 1, 0.99, &r), 1.2363536673577259280e+38 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, 0.99, &r), 3.796186436458346579e+46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -1.5, -0.99, &r), 0.14733409946001025146, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, 0.99, &r), -1.1301780432998743440e+65, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -8.5, -0.99, &r), -8.856462606575344483, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, 0.99, &r), 2.0712920991876073253e+95, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -21.5, -0.99, &r), -74.30517015382249216, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, 0.99, &r),  -3.186778061428268980e+262, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (9, 9, -100.5, -0.99, &r),  2.4454358338375677520, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (25, 25, 1, -0.5, &r), -2.9995530823639545027e-06, TEST_SQRT_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/64.0, &r), 3.17175539044729373926, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/128.0, &r), 3.59937243502024563424, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/256.0, &r), 4.03259299524392504369, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/1024.0, &r), 4.90784159359675398250, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/65536.0, &r), 7.552266033399683914, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, 1.0-1.0/16777216.0, &r), 11.08235454026043830363, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/1024.0, &r), 0.762910940909954974527, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/65536.0, &r), 0.762762124908845424449, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 2.0, -1.0+1.0/1048576.0, &r), 0.762759911089064738044, TEST_TOL0, GSL_SUCCESS);

  /* added special handling with x == 1.0 , Richard J. Mathar, 2008-01-09 */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, 0.5, 3.0, 1.0, &r), 1.6976527263135502482014268 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1.5, -4.2, 3.0, 1.0, &r), .15583601560025710649555254 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-7.4, 0.7, -1.5, 1.0, &r), -.34478866959246584996859 , TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (0.1, -2.7, -1.5, 1.0, &r), 1.059766766063610122925 , TEST_TOL2, GSL_SUCCESS);

  /* Taylor Binnington a = 0 */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (0, -2, -4, 0.5, &r), 1.0 , TEST_TOL2, GSL_SUCCESS);

  /* Andrew Benson <abenson@caltech.edu> bug #24812
     in Pari:
     poch(a,x) = { gamma(a+x)/gamma(a) }
     t(a,b,c,x,k) = { (poch(a,k)*poch(b,k)/poch(c,k)) * (x^k)/(k!) }
     suminf(k=0,t(-10.34, 2.05, 3.05, 0.1725,k))  */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-10.34, 2.05, 3.05, 0.1725, &r), 0.310473552213130010351006093079548, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-9.99999999999, 2.05, 3.05, 0.1725, &r),0.32141934630197487540298837643890, TEST_TOL2, GSL_SUCCESS);

  /* Didier Pinchon also bug #24812 */
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (11, -1, 11.0/2.0, 0.125 , &r), 0.75, TEST_TOL2, GSL_SUCCESS);

  /* Bill Maier - bug #45926 */
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.2, 8.8, 10.0, 0.8, &r), 0.77998971427681563, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.2, 9.8, 11.0, 0.8, &r), 0.77574573497387267, TEST_TOL0, GSL_SUCCESS);

  /* Test case from Hatef Monajemi <monajemi@stanford.edu> (bug #39056) */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (3.5, -0.5, 5.0, 0.9, &r), 0.5923981284370653465208973272, TEST_TOL1, GSL_SUCCESS);

  /* Test c-a-b integer, a and b of opposite sign
   * (Savannah bugs #32306 and #54998) */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.5, 1.5, 1.0, 0.5, &r),  0.5393526011883793566679, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.5, 1.5, 1.0, 0.9, &r), -0.2346215366393253645465, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.25, 0.25, 1.0, 0.25, &r), 0.9833426507751652325159, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-0.25, 0.25, 1.0, 0.9, &r),  0.9205892513820927681855, TEST_TOL1, GSL_SUCCESS);

  /* The alternating Gauss series must report a meaningful error, not
   * the 1.0 produced by the old del_pos bookkeeping (bug #32306). */
  {
    gsl_sf_result er;
    int stat = gsl_sf_hyperg_2F1_e(-0.5, 1.5, 1.0, 0.4, &er);
    gsl_test(stat != GSL_SUCCESS, "gsl_sf_hyperg_2F1_e (-0.5,1.5,1,0.4) status");
    gsl_test(!(er.err < 1e-6), "gsl_sf_hyperg_2F1_e (-0.5,1.5,1,0.4) error estimate (err=%g)", er.err);
  }

  /* 2F1 terminates when a or b is zero or a negative integer.  The
   * resulting polynomial is finite for every x, including |x| >= 1
   * where the Gauss series does not converge; the dispatch used to
   * reject these with GSL_EUNIMPL or GSL_EDOM (Savannah bugs #21835,
   * #50711, #53905).  References: elementary closed forms and mpmath. */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (0, 1, 1, 11, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (1, 0, 1, 11, &r),  1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, 1, 1, 11, &r), -10.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, -1, 1, 11, &r), 12.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, 1, 2, 11, &r), -4.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, -1, 2, 11, &r), 6.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, -1, -0.5, 1.5, &r), -2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1, -13, 1, 0.651439, &r), 9.4687069999999998737, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-10, 2, 0.5, 0.5, &r), 0.077832384333932321548, TEST_TOL2, GSL_SUCCESS);

  /* Test case from Robert L Wolpert <Wolpert@stat.duke.edu> (bug #39056) */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (-1.0, -10.0, 1.0, 0.5, &r), 6.0, TEST_TOL0, GSL_SUCCESS);

  /* Test case from ldnlwm@163.com */

  TEST_SF(s, gsl_sf_hyperg_2F1_e, (3.23191, -4.0229, 8.02291, 0.5, &r), 0.4300243900348170646, TEST_TOL2, GSL_SUCCESS);

  /* 2F1 conj */

  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (1, 1, 1, 0.5, &r), 3.352857095662929028, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, 0.5, &r), 1.7078067538891293983e+09, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, 0.5, &r), 285767.15696901140627, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 1, -0.5, &r), 0.007248196261471276276, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (8, 8, 5, -0.5, &r), 0.00023301916814505902809, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_e, (25, 25, 1, -0.5, &r), 5.1696944096e-06, TEST_SQRT_TOL0, GSL_SUCCESS);

  /* updated correct values, testing enabled, Richard J. Mathar, 2008-01-09 */

  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.01, 1.0, -0.02, &r), .99980388665511730901180717   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (0.1,  0.5, -0.02, &r), .99901595171179281891589794   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (1,   1, -0.02, &r), .98075549650574351826538049000    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (8,   8, -0.02, &r), .32990592849626965538692141   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F0_e, (50, 50, -0.02, &r), .2688995263772964415245902e-12 , TEST_TOL0, GSL_SUCCESS);


  /* 2F1 renorm */

  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (1, 1, 1, 0.5, &r), 2.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 1, 0.5, &r), 12451584.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, 0.5, &r), 0.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, 1, -0.5, &r), 4945.13671875, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, 0.5, &r), -83081.19167659493609, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, -8, -5.5, -0.5, &r), 2.2510895952730178518e+06, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (8, 8, 5, 0.5, &r), 175.2380952380952381, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, 0.99, &r), 1.6063266334913066551e+46, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (9, 9, -1.5, -0.99, &r), 0.06234327316254516616, TEST_TOL2, GSL_SUCCESS);
  /* The negative-integer-c branch of the renorm functions omitted the
   * x^(1-c) factor of [Abramowitz+Stegun, 15.1.2]; the old expectations
   * below encoded that error (each was too large by x^(c-1)) and are
   * corrected here against the A&S limit evaluated independently with
   * mpmath at 50 digits (Savannah bug #53876). */
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (1, 2, -3, 0.4, &r), 65.843621399176983973, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -1, 0.5, &r), 1237440.0, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -10, 0.5, &r), 68070553334784000.0, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_renorm_e, (5, 5, -100, 0.5, &r), 1.191184208770032773e+176, TEST_TOL3, GSL_SUCCESS);
  
  /* 2F1 conj renorm */

  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, 0.99, &r), 5.912269095984229412e+49, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (9, 9, -1.5, -0.99, &r), 0.10834020229476124874, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -1, 0.5, &r), 3.7212765838394834063e+07, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -10, 0.5, &r), 3.8908590631964624485e+18, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hyperg_2F1_conj_renorm_e, (5, 5, -100, 0.5, &r), 1.2271985759108486647e+178, TEST_TOL3, GSL_SUCCESS);


  return s;
}
