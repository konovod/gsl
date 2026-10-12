/* specfunc/airy_osc.h
 *
 * Ai(-t), Bi(-t) and their derivatives for the oscillatory region x < -1
 * are evaluated in double-double.  The old modulus/phase form computed
 * the phase in double precision, whose absolute error grows like
 * eps*t^(3/2) and becomes the relative error of the value (Savannah bug
 * #66808); carrying the phase in double-double removes that growth.
 *
 * For t >= AIRY_OSC_ASYM_MIN the Olver asymptotic expansions
 * (DLMF 9.7.9-9.7.12) are used; their coefficients u_k, v_k come from the
 * exact recurrence 9.7.2.  For 1 < t < AIRY_OSC_ASYM_MIN the Maclaurin
 * series of the Airy equation is summed in double-double, which is
 * accurate there too (the cancellation reaches about 9 digits at t = 10).
 *
 * Included, like cheb_eval_mode.c, into the translation units that need
 * it (airy.c and airy_der.c).
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_trig.h>
#include <math.h>
#include "airy_dd.h"

#define AIRY_OSC_ASYM_MIN 10.0
/* Above this the Dekker split in airy_dd_mul would overflow (it needs
 * |arg| < 2^996); the phase error has long exceeded a radian by then, so
 * the result is reported as a loss instead of being computed. */
#define AIRY_OSC_ASYM_MAX 1.0e150

/* u_k, k = 0..37, from DLMF 9.7.2 (exact rationals rounded to double) */
static double airy_u_data[38] = {
  1.0, 0.06944444444444445, 0.037133487654320986, 0.03799305912780064,
  0.05764919041266972, 0.11609906402551541, 0.2915913992307505, 0.8776669695100169,
  3.079453030173167, 12.341573332345238, 55.62278536591708, 278.46508077760257,
  1533.1694320127956, 9207.206599726414, 59892.51356587907, 419524.87511655106,
  3148257.4178668265, 25198919.871602368, 214288036.96368033, 1929375549.182493,
  18335766937.890568, 183418303528.83255, 1926471158970.4465, 21196999388647.65,
  243826826879716.03, 2926599219297925.0, 3.659030701264313e+16, 4.7576810203630675e+17,
  6.424049357901938e+18, 8.995207427058378e+19, 1.3045132993176097e+21, 1.9570621786581614e+22,
  3.0338710865943386e+23, 4.854832179436167e+24, 8.011464687609593e+25, 1.3621079545263217e+27,
  2.3839516727271057e+28, 4.2915604492858035e+29
};

/* pi/4 and 1/sqrt(pi) to double-double accuracy */
static const airy_dd_t airy_pi_over_4   = { 7.8539816339744828e-01, 3.0616169978683830e-17 };
static const airy_dd_t airy_inv_sqrt_pi = { 5.6418958354775628e-01, 7.6677298065829406e-18 };

/* Ai(0), Ai'(0), Bi(0), Bi'(0) to double-double accuracy */
static const airy_dd_t airy_ai0  = {  3.5502805388781722e-01,  2.0523363243621199e-17 };
static const airy_dd_t airy_ai0p = { -2.5881940379280682e-01,  2.5222431116108321e-17 };
static const airy_dd_t airy_bi0  = {  6.1492662744600068e-01,  5.0899207794891416e-17 };
static const airy_dd_t airy_bi0p = {  4.4828835735382638e-01, -2.5363237774417305e-17 };


/* The auxiliary series of DLMF 9.7.9-9.7.12,
 *
 *   S1 = sum (-1)^k u_{2k}   / z^{2k},   S2 = sum (-1)^k u_{2k+1} / z^{2k+1},
 *   D1 = sum (-1)^k v_{2k}   / z^{2k},   D2 = sum (-1)^k v_{2k+1} / z^{2k+1},
 *
 * with z = zeta and v_k = (6k+1)/(1-6k) u_k.  They are small corrections
 * to cos/sin and do not cancel, so double is enough.  The sums stop at the
 * round-off floor; *serr bounds the first neglected term. */
static int
airy_osc_series(double z, double * S1, double * S2,
                double * D1, double * D2, double * serr)
{
  double inv  = 1.0/z;
  double inv2 = inv*inv;
  double e    = inv2;                    /* z^{-2k}, k >= 1 */
  double s1 = 1.0;
  double s2 = airy_u_data[1]*inv;
  double d1 = 1.0;
  double d2 = (7.0/(1.0 - 6.0))*airy_u_data[1]*inv;   /* v_1 = -7/5 */
  double last = 0.0;
  int k;

  for(k = 1; k <= 18; k++) {
    double sign = (k & 1) ? -1.0 : 1.0;
    double u2k  = airy_u_data[2*k];
    double u2k1 = airy_u_data[2*k+1];
    double v2k  = (12.0*k + 1.0)/(1.0 - 12.0*k)*u2k;
    double v2k1 = (12.0*k + 7.0)/(-12.0*k - 5.0)*u2k1;
    double c1 = sign*u2k*e;
    double c2 = sign*u2k1*e*inv;
    double q1 = sign*v2k*e;
    double q2 = sign*v2k1*e*inv;

    s1 += c1;  s2 += c2;
    d1 += q1;  d2 += q2;
    last = fabs(c1) + fabs(c2);
    e *= inv2;

    if(k >= 2 && last < GSL_DBL_EPSILON)
      break;
  }

  *S1 = s1; *S2 = s2; *D1 = d1; *D2 = d2; *serr = last;
  return GSL_SUCCESS;
}


/* Olver asymptotic expansion of the oscillatory branch, DLMF 9.7.9-9.7.12.
 * deriv = 0 gives Ai(-t)/Bi(-t), deriv = 1 their derivatives; is_bi picks
 * Bi. */
static int
airy_osc_asymp(const double t, const int deriv, const int is_bi,
               double * val, double * err)
{
  airy_dd_t s, z, t14, pref, angle;
  double zeta;
  double S1, S2, D1, D2, serr;
  double ca, sa, bracket, dtheta;
  gsl_sf_result c_hi, s_hi, c_lo, s_lo;
  int stat;

  /* the double-double product cannot represent arguments this large and
   * the phase is indeterminate well before it */
  if(t > AIRY_OSC_ASYM_MAX) {
    t14 = airy_dd_sqrt(airy_dd_sqrt(airy_dd_make(t)));
    *val = 0.0;
    *err = airy_inv_sqrt_pi.hi * (deriv ? t14.hi : 1.0/t14.hi);
    return GSL_ELOSS;
  }

  /* zeta = (2/3) t^{3/2}; angle = zeta - pi/4, both in double-double */
  s     = airy_dd_sqrt(airy_dd_make(t));
  z     = airy_dd_mul(s, airy_dd_make(t));
  z     = airy_dd_div(airy_dd_mul(z, airy_dd_make(2.0)), airy_dd_make(3.0));
  zeta  = z.hi + z.lo;
  angle = airy_dd_sub(z, airy_pi_over_4);

  /* beyond this the phase is not known to a radian and the value is
   * meaningless; report zero with the envelope as the error */
  dtheta = 4.0*GSL_DBL_EPSILON*GSL_DBL_EPSILON*(fabs(zeta) + 1.0);
  if(dtheta >= 1.0) {
    t14 = airy_dd_sqrt(airy_dd_sqrt(airy_dd_make(t)));
    *val = 0.0;
    *err = airy_inv_sqrt_pi.hi * (deriv ? t14.hi : 1.0/t14.hi);
    return GSL_ELOSS;
  }

  airy_osc_series(zeta, &S1, &S2, &D1, &D2, &serr);

  /* cos/sin of angle.hi + angle.lo, the low word carrying the residual */
  stat  = gsl_sf_cos_e(angle.hi, &c_hi);
  stat |= gsl_sf_sin_e(angle.hi, &s_hi);
  stat |= gsl_sf_cos_e(angle.lo, &c_lo);
  stat |= gsl_sf_sin_e(angle.lo, &s_lo);
  ca = c_hi.val*c_lo.val - s_hi.val*s_lo.val;
  sa = s_hi.val*c_lo.val + c_hi.val*s_lo.val;

  if(deriv)
    bracket = is_bi ? ( ca*D1 + sa*D2) : ( sa*D1 - ca*D2);
  else
    bracket = is_bi ? (-sa*S1 + ca*S2) : ( ca*S1 + sa*S2);

  /* prefactor t^{1/4}/sqrt(pi) (derivative) or its reciprocal (value),
   * in double-double so it does not limit the accuracy */
  t14 = airy_dd_sqrt(airy_dd_sqrt(airy_dd_make(t)));
  pref = deriv ? airy_dd_mul(airy_inv_sqrt_pi, t14)
               : airy_dd_div(airy_inv_sqrt_pi, t14);

  {
    airy_dd_t r = airy_dd_mul(airy_dd_make(bracket), pref);
    *val = r.hi + r.lo;
  }

  *err  = fabs(pref.hi) * (dtheta + serr
                           + fabs(c_hi.err) + fabs(s_hi.err)
                           + fabs(c_lo.err) + fabs(s_lo.err));
  *err += 2.0 * GSL_DBL_EPSILON * fabs(*val);
  if(*err == 0.0) *err = GSL_DBL_EPSILON;

  return stat;
}


/* Maclaurin series of the Airy equation in double-double,
 *
 *   T[0] = Ai(0) etc.,  T[n] = T[n-3]/(n(n-1)),
 *   Ai(x) = sum T[n] x^n,   Ai'(x) = sum n T[n] x^{n-1},
 *
 * used for 1 < t < AIRY_OSC_ASYM_MIN, where the cancellation is small
 * enough for double-double. */
static int
airy_osc_maclaurin(const double t, const int deriv, const int is_bi,
                   double * val, double * err)
{
  const double x  = -t;
  const airy_dd_t c0 = is_bi ? airy_bi0  : airy_ai0;
  const airy_dd_t c1 = is_bi ? airy_bi0p : airy_ai0p;
  const int nmax = 200;
  airy_dd_t T[203];
  airy_dd_t xpow, xd, sum, dsum;
  double maxterm = 0.0;
  int n, last = 0;
  int small = 0;

  T[0] = c0;
  T[1] = c1;
  T[2] = airy_dd_make(0.0);
  for(n = 3; n <= nmax; n++)
    T[n] = airy_dd_div(T[n-3], airy_dd_make((double)n*(double)(n-1)));

  xd   = airy_dd_make(x);
  xpow = airy_dd_make(1.0);
  sum  = airy_dd_make(0.0);
  dsum = airy_dd_make(0.0);

  for(n = 0; n <= nmax; n++) {
    airy_dd_t term = airy_dd_mul(T[n], xpow);
    double mag = fabs(term.hi) + fabs(term.lo);

    if(n >= 1) {
      airy_dd_t dterm = airy_dd_mul(airy_dd_make((double)n), term);
      double m;
      dterm = airy_dd_div(dterm, xd);       /* n T[n] x^{n-1} */
      dsum = airy_dd_add(dsum, dterm);
      m = fabs(dterm.hi) + fabs(dterm.lo);
      if(m > mag) mag = m;
    }
    sum = airy_dd_add(sum, term);
    if(mag > maxterm) maxterm = mag;
    last = n;

    /* The recurrence makes every third T[n] vanish, so a single small term
     * is not evidence of convergence; wait until three in a row are at the
     * round-off floor. */
    if(n >= 30 && mag <= GSL_DBL_EPSILON*GSL_DBL_EPSILON*maxterm)
      small++;
    else
      small = 0;
    if(small >= 3)
      break;

    xpow = airy_dd_mul(xpow, xd);
  }

  {
    airy_dd_t r = deriv ? dsum : sum;
    *val = r.hi + r.lo;
  }

  *err  = 4.0*GSL_DBL_EPSILON*GSL_DBL_EPSILON*maxterm*(double)(last+1);
  *err += 2.0 * GSL_DBL_EPSILON * fabs(*val);
  if(*err == 0.0) *err = GSL_DBL_EPSILON;

  return GSL_SUCCESS;
}
