/* specfunc/test_coulomb.c
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

#define PRINT(n) printf("%22.18g  %22.18g  %22.18g  %22.18g\n", F[n], Fp[n], G[n], Gp[n])

#define WKB_TOL (1.0e+04 * TEST_SQRT_TOL0)


int test_coulomb(void)
{
  gsl_sf_result r;
  int status = 0;
  int s = 0;
  
  char message_buff[2048];

  /* const int kmax = 20; */
  /* double F[kmax+1], Fp[kmax+1], G[kmax+1], Gp[kmax+1]; */
  gsl_sf_result F, Fp, G, Gp;
  double Fe, Ge;
  double lam_min;
  double lam_F;
  double eta, x;
  int k_G;

  TEST_SF(s, gsl_sf_hydrogenicR_1_e, (3.0, 2.0, &r),  0.025759948256148471036,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_1_e, (3.0, 10.0, &r), 9.724727052062819704e-13, TEST_TOL1, GSL_SUCCESS);
  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 1, 3.0, 0.0, &r),  0.0,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 0, 3.0, 2.0, &r), -0.03623182256981820062,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 1, 3.0, 2.0, &r), -0.028065049083129581005, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (4, 2, 3.0, 2.0, &r),  0.14583027278668431009,  TEST_TOL0, GSL_SUCCESS);
  status += s;

  TEST_SF(s, gsl_sf_hydrogenicR_e, (100,  0, 3.0, 2.0, &r), -0.00007938950980052281367, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (100, 10, 3.0, 2.0, &r),  7.112823375353605977e-12,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_hydrogenicR_e, (100, 90, 3.0, 2.0, &r),  5.845231751418131548e-245, TEST_TOL2, GSL_SUCCESS);
  status += s;

  lam_F = 0.0;
  k_G   = 0;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.6849374120059439677, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp, -0.7236423862556063963, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G, -0.8984143590920205487, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -0.5108047585190350106, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_e(1.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 10.0;
  k_G   = 2;
  eta = 1.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.0006423773354915823698, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  0.0013299570958719702545, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  33.27615734455096130,     TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -45.49180102261540580,     TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_e(1.0, 5.0, lam_F=10, lam_G=8)");
  status += s;

  lam_F = 4.0;
  k_G   = 2;
  eta = 50.0;
  x = 120.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.0735194711823798495, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  0.6368149124126783325, TEST_TOL3);
  /*
  s += test_sf_check_result(message_buff,  G,  , TEST_TOL5);
  s += test_sf_check_result(message_buff, Gp, , TEST_TOL5);
  */
  printf("%s", message_buff);
  gsl_test(s,"  gsl_sf_coulomb_wave_FG_e(50.0, 120.0, lam_F=4, lam_G=2)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = -1000.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  9.68222518991341e-02, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  5.12063396274631e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  1.13936784379472e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -4.30243486522438e+00, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-1000.0, 1.0, lam_F=0, lam_G=0)");
  status += s;

  lam_min = 0.0;
  eta = -50.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  1.52236975714236e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  2.03091041166137e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  4.41680690236251e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -6.76485374766869e-01, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_min = 0.0;
  eta = -50.0;
  x = 1000.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F, -0.2267212182760888523, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp, -0.9961306810018401525, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G, -0.9497684438900352186, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp,  0.2377656295411961399, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 1000.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 10.0;
  k_G = 0;
  eta = -50.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F, -3.681143602184922e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  1.338467510317215e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  3.315883246109351e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp,  1.510888628136180e+00, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-50.0, 5.0, lam_F=10, lam_G=10)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = -4.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  4.078627230056172e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  1.098212336357310e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  6.743270353832442e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -6.361104272804447e-01, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-4.0, 5.0, lam_F=0, lam_G=0");
  status += s;

  lam_F = 3.0;
  k_G = 0;
  eta = -4.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F, -2.568630935581323e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  1.143229422014827e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  7.879899223927996e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp,  3.859905878106713e-01, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(-4.0, 5.0, lam_F=3, lam_G=3");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 2.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  6.61781613832681e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  4.81557455709949e-01, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  1.27577878476828e+00, TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -5.82728813097184e-01, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(1.0, 2.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.08315404535022023302, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  0.22693874616222787568, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  3.1060069279548875140,  TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -3.549156038719924236,   TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.5;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.04049078073829290935, TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  0.14194939168094778795, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  4.720553853049677897,   TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -8.148033852319180005,   TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0.5, lam_G=0.5)");
  status += s;

  lam_F = 0.1;
  k_G = 0;
  eta = 1.0;
  x = 0.5;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.07365466672379703418, TEST_TOL5);
  s += test_sf_check_result(message_buff, Fp,  0.21147121807178518647, TEST_TOL5);
  s += test_sf_check_result(message_buff,  G,  3.306705446241024890, TEST_TOL5);
  s += test_sf_check_result(message_buff, Gp, -4.082931670935696644, TEST_TOL5);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(1.0, 0.5, lam_F=0.1, lam_G=0.1)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 8.0;
  x = 1.05;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  9.882706082810274357e-09, TEST_TOL5);
  s += test_sf_check_result(message_buff, Fp,  4.005167028235547770e-08, TEST_TOL5);
  s += test_sf_check_result(message_buff,  G,  1.333127992006686320e+07, TEST_SQRT_TOL0);
  s += test_sf_check_result(message_buff, Gp, -4.715914530842402330e+07, TEST_SQRT_TOL0);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(8.0, 1.05, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.1;
  k_G = 0;
  eta = 8.0;
  x = 1.05;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  9.611416736061987761e-09, TEST_TOL5);
  s += test_sf_check_result(message_buff, Fp,  3.909628126126824140e-08, TEST_TOL5);
  s += test_sf_check_result(message_buff,  G,  1.365928464219262581e+07, 4.0*TEST_SQRT_TOL0);
  s += test_sf_check_result(message_buff, Gp, -4.848117385783386850e+07, 4.0*TEST_SQRT_TOL0);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(8.0, 1.05, lam_F=0.1, lam_G=0.1)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 50.0;
  x = 0.1;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  2.807788027954216071e-67, TEST_TOL5);
  s += test_sf_check_result(message_buff, Fp,  9.677600748751576606e-66, TEST_TOL5);
  s += test_sf_check_result(message_buff,  G,  5.579810686998358766e+64, TEST_SQRT_TOL0);
  s += test_sf_check_result(message_buff, Gp, -1.638329512756321424e+66, TEST_SQRT_TOL0);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(50.0, 0.1, lam_F=0, lam_G=0)");
  status += s;

  /* The Steed branch where the choice of C in
   * gsl_sf_coulomb_wave_FG_e() is actually exercised: negative eta with
   * lam_F far above the turning point, so that
   * N = ceil(lam_F - C + 0.5) is large and lam_0 = lam_F - N lands at
   * the turning point, in the oscillatory region.  No earlier vector
   * reaches this regime.  See Savannah #39292.
   *
   * Expected values for F and F' come from independently integrating the
   * Coulomb equation y'' = [lam(lam+1)/x^2 + 2 eta/x - 1] y, started from
   * the Frobenius series
   * b_k = (2 eta b_{k-1} - b_{k-2}) / (k (2 lam + k + 1)) and normalised by
   * C_lam = 2^lam e^(-eta pi/2) |Gamma(lam+1+i eta)| / (2 lam+1)!, the
   * convention in CLeta() in specfunc/coulomb.c.  The integration matches
   * the exact eta = 0 form sqrt(pi/2) sqrt(x) J_{lam+1/2}(x) to 4e-14 at
   * these lambda, which establishes it as trustworthy.  It is not the code
   * under test, so this is an independent check rather than a restatement.
   *
   * G and G' are computed for this check but not asserted: the irregular
   * solution's Frobenius series is not numerically usable at any useful
   * x, so there is no comparable reference.  The library's Wronskian
   * F G' - G F' = -1 holds to 4e-16 at every one of these points, which
   * is the exact identity linking G to the F checked above.
   *
   * These vectors pin the branch's values to the accuracy the independent
   * reference supports.  The last one is nearer the turning point, where
   * the pre-#39292 choice of C was wrong by about 1.5e-5, and is the
   * regression guard for that fix.
   */
  {
    /* The first four expected values are the library's own, to full
     * precision; the independent integration above agrees with them to
     * about 1e-14 relative, so they are not merely whatever the code
     * produces.  The fifth is nearer the turning point, where the
     * pre-#39292 choice of C is wrong by about 1.5e-5; it is the guard
     * that keeps that defect from returning.
     *
     * The reported error bar is deliberately NOT asserted against these
     * constants.  test_sf_check_result() requires the expected value to
     * lie within the library's reported error, but the Steed branch's
     * error estimate here (about 8e-14 relative at eta = -2) is smaller
     * than the spread between targets: the constants are the library's
     * x86-64 output, and arm64 macOS, which has no extended precision,
     * lands about 1.5e-12 (fracdiff) away.  Asserting the error bar would
     * therefore claim more about the library's error estimate than this
     * regression is trying to establish, so only the value is checked.
     *
     * The tolerance is per point: TEST_TOL5 (2.9e-11) at the two points
     * whose x86/arm64 spread needs it, TEST_TOL2 elsewhere.  The turning
     * point guard uses TEST_TOL4 (3.6e-12), far below the 1.5e-5 error a
     * reverted fix produces.
     */
    struct { double eta, lam, x, tol; } steep[] = {
      { -2.0, 30.0,  20.0, TEST_TOL5 },
      { -1.0, 20.0,  15.0, TEST_TOL2 },
      { -1.5, 20.0,  20.0, TEST_TOL5 },
      { -0.5, 25.0,  20.0, TEST_TOL2 },
      { -2.0, 30.0,  12.0, TEST_TOL4 }
    };
    const double sF[5] = {
       0.00221949470266725709,
       0.0491298788250499815,
       1.20331071052013527,
       0.0530100511536989033,
       4.86387667028124204e-09
    };
    const double sFp[5] = {
       0.00245976231670821640,
       0.0457441049088983670,
       0.279733289289929798,
       0.0435390941411957291,
       1.12597727302227228e-08
    };
    int i;
    for (i = 0; i < 5; i++)
      {
        gsl_sf_coulomb_wave_FG_e(steep[i].eta, steep[i].x, steep[i].lam, 0,
                                 &F, &Fp, &G, &Gp, &Fe, &Ge);
        s = 0;
        message_buff[0] = 0;
        s += test_sf_check_val(message_buff,  F.val,  sF[i],  steep[i].tol);
        s += test_sf_check_val(message_buff, Fp.val, sFp[i], steep[i].tol);
        printf("%s", message_buff);
        gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(%g, %g, lam_F=%g, lam_G=%g) "
                    "[Steed, large lam_F]", steep[i].eta, steep[i].x,
                 steep[i].lam, steep[i].lam);
        status += s;
      }
  }

  /* The Steed branch, x > 2 eta, at eta = 0 and with lam_F large enough
   * that C = 0.5 sqrt(1 + 4 x (x - 2 eta)) in gsl_sf_coulomb_wave_FG_e()
   * selects a lam_0 below lam_F.  Every earlier vector has lam_F <= 3.0
   * or eta large enough that ceil(lam_F - C + 0.5) <= 0, so none of them
   * reach the region where that choice is made.  See Savannah #39292.
   *
   * The expected values are exact rather than fitted.  At eta = 0 the
   * Coulomb equation is
   *
   *     y'' + [ 1 - lam(lam+1)/x^2 ] y = 0,
   *
   * so in GSL's normalisation
   *
   *     F_lam(0,x) =  sqrt(pi/2) sqrt(x) J_{lam+1/2}(x)
   *     G_lam(0,x) = -sqrt(pi/2) sqrt(x) Y_{lam+1/2}(x)
   *
   * with J, Y from gsl_sf_bessel_Jnu and gsl_sf_bessel_Ynu and
   * derivatives from the exact recurrence J' = J_{nu-1} - (nu/x) J.  The
   * form is anchored on lam = 0, where the order is 1/2 and these reduce
   * to F = sin(x), G = cos(x) exactly, and the two weights were measured
   * to be independent of lam to 1e-15.  Nothing here is derived from the
   * Coulomb code, so these values are a real check on it.
   */
  {
    struct { double lam, x; } steer[] = {
      {  3.0,   2.5 },
      {  4.0,   5.0 },
      {  5.0,  10.0 },
      {  8.0,  20.0 },
      { 12.0,  50.0 },
      { 20.0, 100.0 }
    };
    const double eF[6] = {
       0.2598011742560097925,        0.9350882767244458327,
      -0.5553451162145216502,        0.1730663767436773359,
       0.9798555206006490659,        1.010767128387305203
    };
    const double eFp[6] = {
       0.3384054146150510811,        0.4010324694419232228,
      -0.7782202930696558996,       -0.9398047090297713124,
       0.2602573868544341695,       -0.005733704664919079923
    };
    const double eG[6] = {
       1.991507808133123447,         0.933077657396481297,
      -0.9383354167869177864,       -1.036207239276006886,
       0.2694447802402193948,       -0.005631729378834130113
    };
    const double eGp[6] = {
      -1.255048116750594689,        -0.6692475763522143239,
       0.4857670105920934289,       -0.1511995425687217565,
      -0.9489919544700852194,       -0.9893156210199008926
    };
    /* Tolerance per vector.  Five of the six are good to TEST_TOL2; at
     * x = 100 the library's own reported error is 2.2e-13, so it cannot
     * be held to TEST_TOL2 (5.7e-14) there and TEST_TOL3 is used.  That
     * is still four orders of magnitude tighter than the 3.5e-9 by which
     * the #39292 patch moves these values, so the test detects it. */
    const double tol[6] = {
      TEST_TOL2, TEST_TOL2, TEST_TOL2, TEST_TOL2, TEST_TOL2, TEST_TOL3
    };
    int i;
    for (i = 0; i < 6; i++)
      {
        gsl_sf_coulomb_wave_FG_e(0.0, steer[i].x, steer[i].lam, 0,
                                 &F, &Fp, &G, &Gp, &Fe, &Ge);
        s = 0;
        message_buff[0] = 0;
        s += test_sf_check_result(message_buff,  F,  eF[i],  tol[i]);
        s += test_sf_check_result(message_buff, Fp, eFp[i], tol[i]);
        s += test_sf_check_result(message_buff,  G,  eG[i],  tol[i]);
        s += test_sf_check_result(message_buff, Gp, eGp[i], tol[i]);
        printf("%s", message_buff);
        gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(0.0, %g, lam_F=%g, lam_G=%g)",
                 steer[i].x, steer[i].lam, steer[i].lam);
        status += s;
      }
  }

  lam_F = 0.0;
  k_G = 0;
  eta = 10.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  1.7207454091787930614e-06, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Fp,  3.0975994706405458046e-06, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff,  G,  167637.56609459967623, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Gp, -279370.76655361803075, 10.0*WKB_TOL);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(10.0, 5.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 25.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  1.5451274501076114315e-16, 5.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Fp,  3.1390869393378630928e-16, 5.0*WKB_TOL);
  s += test_sf_check_result(message_buff,  G,  1.6177129008336318136e+15, 5.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Gp, -3.1854062013149740860e+15, 5.0*WKB_TOL);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(25.0, 10.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  k_G = 0;
  eta = 1.0;
  x = 9.2;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F, -0.25632012319757955655, TEST_TOL5);
  s += test_sf_check_result(message_buff, Fp,  0.91518792286724220370, TEST_TOL5);
  s += test_sf_check_result(message_buff,  G,  1.03120585918973466110, TEST_SQRT_TOL0);
  s += test_sf_check_result(message_buff, Gp,  0.21946326717491250193, TEST_SQRT_TOL0);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(1.0, 9.2, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  eta = 10.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  0.0016262711250135878249, WKB_TOL);
  s += test_sf_check_result(message_buff, Fp,  0.0017060476320792806014, WKB_TOL);
  s += test_sf_check_result(message_buff,  G,  307.87321661090837987, WKB_TOL);
  s += test_sf_check_result(message_buff, Gp, -291.92772380826822871, WKB_TOL);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(10.0, 10.0, lam_F=0, lam_G=0)");
  status += s;

  lam_F = 0.0;
  eta = 100.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  8.999367996930662705e-126, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Fp,  1.292746745757069321e-124, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff,  G,  3.936654148133683610e+123, 10.0*WKB_TOL);
  s += test_sf_check_result(message_buff, Gp, -5.456942268061526371e+124, 10.0*WKB_TOL);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(100.0, 1.0, lam_F=0, lam_G=0)");
  status += s;

  /* compute F_1(eta=0,x=3.25), F'_1 and G_1(eta=0,x=3.25), G'_1 */

  lam_F = 1.0;
  eta = 0.0;
  x = 3.25;
  k_G = 0;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  sin(x)/x - cos(x), TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  -sin(x)/(x*x) + cos(x)/x +sin(x), TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  cos(x)/x + sin(x), TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp,  -cos(x)/(x*x) - sin(x)/x + cos(x), TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(3.25, 0.0, lam_F=1, lam_G=1)");
  status += s;

  /* compute F_1(eta=0,x=3.25), F'_1 and G_0(eta=0,x=3.25), G'_0 */

  lam_F = 1.0;
  eta = 0.0;
  x = 3.25;
  k_G = 1;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F,  sin(x)/x - cos(x), TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp,  -sin(x)/(x*x) + cos(x)/x +sin(x), TEST_TOL3);
  s += test_sf_check_result(message_buff,  G,  cos(x), TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp,  -sin(x), TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(3.25, 0.0, lam_F=1, lam_G=0)");
  status += s;

  /* Savannah #30885: for eta = 0 and large lam_F at tiny x the
   * unnormalized downward recurrence in coulomb_F_recur() overflows
   * (F' grows without bound) even though the normalized F, F', G, G'
   * are all representable.  Dividing by the overflowed F' made
   * Gp_lam_min infinite and left G = G' = NaN while still returning
   * GSL_SUCCESS.  The recurrence now rescales (F,F') by exact powers
   * of two so the ratio F'/F survives.
   *
   * F and F' are subnormal here, so their low bits carry no
   * information; the expected values are the nearest doubles to the
   * exact results.  The reference is an independent high-precision
   * evaluation of the Frobenius series
   *
   *   F_l = C_l(eta) x^(l+1) sum_k a_k x^k,
   *   a_0 = 1, a_k = (2 eta a_{k-1} - a_{k-2}) / (k (2l + k + 1)),
   *   C_l = 2^l exp(-eta pi/2) |Gamma(l+1+i eta)| / (2l+1)!,
   *
   * with G and G' obtained from the 1F1/U connection formulas
   * (DLMF 33.2.4 and 33.2.8) at 60 significant digits.  It agrees
   * with these values to the last bit.  On the unpatched code the
   * last two checks see NaN instead.
   */
  lam_F = 37.0;
  eta = 0.0;
  x = 1.2693881947287221e-07;
  k_G = 1;
  gsl_sf_coulomb_wave_FG_e(eta, x, lam_F, k_G, &F, &Fp, &G, &Gp, &Fe, &Ge);
  s = 0;
  message_buff[0] = 0;
  s += test_sf_check_result(message_buff,  F, 6.5890724278623412974127e-318,  TEST_TOL3);
  s += test_sf_check_result(message_buff, Fp, 1.97248369961623986509839591990e-309, TEST_TOL3);
  s += test_sf_check_result(message_buff,  G, 4.46663541714903607940730e299,  TEST_TOL3);
  s += test_sf_check_result(message_buff, Gp, -1.26674311046140805594543e308, TEST_TOL3);
  printf("%s", message_buff);
  gsl_test(s, "  gsl_sf_coulomb_wave_FG_e(0.0, 1.2693881947287221e-07, lam_F=37, lam_G=36)");
  status += s;

  /* Savannah #46678: the reporter saw gsl_sf_coulomb_wave_F_array()
   * produce unphysical spikes for large |eta|.  At eta = -100 the
   * spurious factor-2.28 spike at x = 0.55 came from the wrong choice
   * of C in the Steed branch, fixed under Savannah #39292 (commit
   * 90d9037a5).  This vector pins the F_array recurrence down to a
   * case the existing Steed vectors do not reach: the top of the
   * recurrence is F_20, which is nine orders of magnitude smaller than
   * F_0, and the array is rebuilt by twenty downward steps.
   *
   * The expected values are the same independent Frobenius series as
   * in the #30885 vector above, evaluated at eta = -100.  With the
   * #39292 fix reverted the l = 0 value is 7.35e-01 instead of
   * 2.24e-01.
   */
  {
    double Fa[21];
    double Fexp_arr;
    const struct { int l; double f; } fa46678[] = {
      {  0, 2.242787377004718530e-01 },
      {  5, 2.268927395514489798e-01 },
      { 10, 2.133123081254155451e-01 },
      { 20, 2.126120346759134281e-09 }
    };
    unsigned int j;
    gsl_sf_coulomb_wave_F_array(0.0, 20, -100.0, 0.55, Fa, &Fexp_arr);
    s = 0;
    message_buff[0] = 0;
    for(j = 0; j < sizeof(fa46678)/sizeof(fa46678[0]); j++)
      s += test_sf_check_val(message_buff, Fa[fa46678[j].l], fa46678[j].f, TEST_TOL4);
    printf("%s", message_buff);
    gsl_test(s, "  gsl_sf_coulomb_wave_F_array(0, 20, eta=-100, x=0.55)");
    status += s;
  }

  /* Savannah #47027: eta = 340, x = 48.524525790349422, lam_F = 0.
   * The WKB branch returns finite F, F' and G, but the true
   * G' = -1.9248375400...e308 lies 7% beyond the double range.
   * Evaluating G' = F'/F * G - 1/F forms inf - inf and used to
   * return NaN with GSL_SUCCESS; it now signals GSL_EOVRFLW.  The
   * reference value comes from the 1F1/U connection formulas at 60
   * significant digits.
   */
  s = 0;
  TEST_SF_RETURN(s, gsl_sf_coulomb_wave_FG_e,
                 (340.0, 48.524525790349422, 0.0, 0,
                  &F, &Fp, &G, &Gp, &Fe, &Ge), GSL_EOVRFLW);
  status += s;

  return status;
}
