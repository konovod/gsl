/* specfunc/airy_dd.h
 *
 * Minimal double-double arithmetic for the Airy oscillatory branch.
 *
 * The modulus/phase representation of Ai(-t) and Bi(-t) computes the
 * phase in double precision, so its absolute error grows like
 * eps*t^(3/2); that error is the relative error of the value.  Carrying
 * the phase (and the series that build it) in double-double removes the
 * growth.  Two doubles hold about 31 significant digits.
 *
 * The product is Dekker's exact split product (Dekker 1971).  It needs
 * no FMA and is correct for any rounding mode; a fused multiply-add only
 * makes the residual exact, which is harmless here.
 *
 * These helpers are `static inline` and included, like cheb_eval_mode.c,
 * into the translation units that need them.
 */

#ifndef GSL_SPECFUNC_AIRY_DD_H
#define GSL_SPECFUNC_AIRY_DD_H

#include <math.h>

typedef struct { double hi, lo; } airy_dd_t;

static inline airy_dd_t
airy_dd_make(double x)
{
  airy_dd_t r;
  r.hi = x;
  r.lo = 0.0;
  return r;
}

/* s + e == a.hi + b.hi with |e| <= 0.5 ulp(s); the branch picks the
 * larger addend so that the Fast2Sum error term is exact. */
static inline airy_dd_t
airy_dd_add(airy_dd_t a, airy_dd_t b)
{
  airy_dd_t r;
  double s, e;
  s = a.hi + b.hi;
  if(fabs(a.hi) > fabs(b.hi))
    e = (a.hi - s) + b.hi;
  else
    e = (b.hi - s) + a.hi;
  e += a.lo + b.lo;
  r.hi = s + e;
  r.lo = e - (r.hi - s);
  return r;
}

static inline airy_dd_t
airy_dd_neg(airy_dd_t a)
{
  airy_dd_t r;
  r.hi = -a.hi;
  r.lo = -a.lo;
  return r;
}

static inline airy_dd_t
airy_dd_sub(airy_dd_t a, airy_dd_t b)
{
  return airy_dd_add(a, airy_dd_neg(b));
}

/* a*b == p + e exactly, via the Dekker split (2^27 + 1 splits a double
 * into two 26-bit halves whose products are exact). */
static inline void
airy_dd_two_prod(double a, double b, double * p, double * e)
{
  const double split = 134217729.0;  /* 2^27 + 1 */
  double ca, ahi, alo, cb, bhi, blo;

  ca  = split * a;
  ahi = ca - (ca - a);
  alo = a - ahi;
  cb  = split * b;
  bhi = cb - (cb - b);
  blo = b - bhi;

  *p = a * b;
  *e = ((ahi*bhi - *p) + ahi*blo + alo*bhi) + alo*blo;
}

static inline airy_dd_t
airy_dd_mul(airy_dd_t a, airy_dd_t b)
{
  airy_dd_t r;
  double p, e, t;

  airy_dd_two_prod(a.hi, b.hi, &p, &e);
  e += a.hi*b.lo + a.lo*b.hi;
  t = p + e;
  r.hi = t;
  r.lo = e - (t - p);
  return r;
}

static inline airy_dd_t
airy_dd_div(airy_dd_t a, airy_dd_t b)
{
  airy_dd_t q, p, rem, negp;
  double q1, q2;
  airy_dd_t r;

  q1 = a.hi / b.hi;
  q = airy_dd_make(q1);
  p = airy_dd_mul(q, b);
  negp.hi = -p.hi;
  negp.lo = -p.lo;
  rem = airy_dd_add(a, negp);
  q2 = (rem.hi + rem.lo) / b.hi;
  r.hi = q1 + q2;
  r.lo = (q1 - r.hi) + q2;
  return r;
}

/* sqrt(a) in double-double; one Newton step from the double sqrt, whose
 * relative error eps is squared by the correction. */
static inline airy_dd_t
airy_dd_sqrt(airy_dd_t a)
{
  double s = sqrt(a.hi);
  airy_dd_t p, rem;
  double corr;

  p = airy_dd_mul(airy_dd_make(s), airy_dd_make(s));
  rem = airy_dd_sub(a, p);
  corr = (rem.hi + rem.lo) / (2.0 * s);
  return airy_dd_add(airy_dd_make(s), airy_dd_make(corr));
}

#endif /* GSL_SPECFUNC_AIRY_DD_H */
