# Savannah open-issue triage

Generated 2026-10-02 from the `scripts/savannah_bugs.py` sweep (`inventory.tsv`).

Scratch index, not part of the fork's record of changes. Verdicts are curated from
`FORKNEWS`, `SAVANNAH_REVIEW.md` and the git history; every other row is un-triaged.

**219 open items: 123 reviewed/handled, 96 remaining.**

| status | count |
|---|---:|
| fixed | 83 |
| partial | 3 |
| rejected | 33 |
| deferred | 3 |
| superseded | 1 |
| not reviewed (feature-shaped) | 21 |
| not reviewed | 75 |

Legend for the *Patch* column: the sweep's `git apply` verdict (`clean` / `partial` / `dirty`), `source` for a whole
replacement file, `inline` for a patch pasted into the bug text, `none` for no patch.

| # | Date | Kind | Cat | Patch | Fork status | Notes |
|---|---|---|---|---|---|---|
| 21828 | 2007-12-18 | - | Performance | none | not reviewed | suboptimal performance of gsl_fdfsolver_lmsder |
| 21831 | 2007-12-18 | - | Accuracy problem | source levy.c | not reviewed | Levý random number generator for alpha < 1 |
| 21833 | 2007-12-18 | - | Performance | none | not reviewed | suboptimal performance of gsl permutation? |
| 21835 | 2007-12-18 | - | Accuracy problem | clean test_hyperg.diff +1 | **partial** — non-positive-integer termination for x >= 1 fixed (2aa89bae8); c = a+b near x = 1 still returns GSL_EMAXITER | gsl_sf_hyperg_2F1 problematic arguments |
| 21836 | 2007-12-18 | - | Accuracy problem | none | **fixed** — Q returns exact complement of P in the series window (6f23a4cb5) | gamma_inc_P and gamma_inc_Q only satisfy P+Q=1 within errors |
| 21837 | 2007-12-18 | - | Runtime error | none | not reviewed | gsl_linalg_solve_symm_tridiag requires positive definite matrix |
| 24252 | 2008-09-12 | feature | - | source gamma_tail_jpl_080908.c | not reviewed | suggestion: add gamma tail distribution |
| 24871 | 2008-11-18 | feature | - | none | not reviewed | suggestion, add support for E_n |
| 25320 | 2009-01-14 | - | Accuracy problem | none | not reviewed | Import fresnel, bugs on GSL Extension Fresnel |
| 28267 | 2009-12-11 | - | Accuracy problem | source hyperg1F1.c | **partial** — same defect as #43809 (52505315d); transition region `x ~ a^2` still loses digits | poor convergence region for gsl_sf_hyperg_1F1 |
| 29834 | 2010-05-09 | - | Runtime error | source error_cblas_v2.h | not reviewed | insufficient argument checking in blas wrapper |
| 30324 | 2010-07-02 | - | Accuracy problem | none | **deferred** — feature: extend 2F1 to x < -1 by transformation; new domain/algorithm, out of eligibility | improve range of 2F1 |
| 30510 | 2010-07-21 | - | Runtime error | none | **fixed** — U(a,b,x) for x < 0, integer b, non-integer a, via the DLMF 13.2.9 limit (00859f816) | problems with hyperg_U(a,b,x) for x<0 |
| 30540 | 2010-07-24 | feature | Accuracy problem | partial bug-ode2.c +2 | not reviewed | please, add convergence checks in rk4imp/rk2imp |
| 30583 | 2010-07-28 | doc | Documentation | none | **fixed** — Legendre/Carlson relations and the negative-parameter (imaginary-modulus) transformation documented (b0eec8bc6) | improve documentation for Elliptic functions |
| 30885 | 2010-08-27 | - | Runtime error | none | **fixed** — Coulomb F recurrence rescaled, no overflow (25841970e) | nans from gsl_sf_coulomb_wave_FG_e(1.2693881947287221e-07, 0.0, lam_F=37, lam_G=36) |
| 30947 | 2010-09-02 | - | - | clean 0001-Fixed-step-size-control-object.patch | not reviewed | Please, include fixed step size control object for ode suite |
| 31109 | 2010-09-23 | - | Performance | none | not reviewed | ode-initval/bsimp is always high order |
| 31362 | 2010-10-18 | bug | Runtime error | none | **fixed** — NaN rejected as GSL_EDOM in the complete elliptic integrals (14d595eb8) | The Complete Elliptic Integrals (gsl_sf_ellint_Ecomp and _Kcomp) Loop Forever with NaN Argument |
| 31426 | 2010-10-23 | - | Runtime error | none | not reviewed | infinite loop in gsl_eigen_symm |
| 32257 | 2011-01-26 | - | - | none | not reviewed | RFE: Import integration routines from quadrule |
| 32306 | 2011-01-31 | bug | Accuracy problem | source hyp.c | **fixed** — integer-d 2F1 series and error estimate (e4c4ac326, 882c8361d) | sign error in gsl_sf_hyperg_2F1 |
| 32776 | 2011-03-14 | feature | - | source quadratic.c | **fixed** — multimin quadratic minimiser (53cd0d098) | RFE: Add brute-force quadratic numerical multidimensional minimizer |
| 34361 | 2011-09-22 | - | Runtime error | none | not reviewed | gsl_bspline_knots_greville needs inequality constrained linear least squares |
| 35032 | 2011-12-11 | doc | Documentation | none | **fixed** — gsl_test support functions documented in the usage chapter, marked as an unreviewed AI draft (8ef88af6e) | gsl_test.h lacks documentation in the reference manual |
| 36152 | 2012-04-11 | bug | Accuracy problem | source testbessel.c | **fixed** — Y family libm sin/cos (a08ef2f7f) + exact sin/cos reduction (8a46ec7cf) | Incorrect asymptotics of spherical Bessel functions |
| 36197 | 2012-04-15 | - | Build | dirty 36197b.diff +1 | not reviewed | reserved identifier violation |
| 36578 | 2012-06-02 | - | Documentation | none | **fixed** — doc now says the error handler is invoked with GSL_EUNSUP and that the void function returns no code (060c2472c) | inconsistency in gsl_ieee_env_setup doc/api |
| 37209 | 2012-08-28 | bug | Accuracy problem | none | **fixed** — upstream 441bc40ff (inherited); NaN gone | gsl_sf_bessel_jl_e returns NaN for large inputs without setting error code |
| 37408 | 2012-09-20 | - | - | clean au.patch | **fixed** — const AU/parsec updated to IAU 2012 (3377c7fbe) | The astronomical unit (AU) has been re-defined |
| 37894 | 2012-12-10 | bug | Build | clean gsl-autotools.diff | **fixed** — upstream 1d002ee93/ae19e3e8b (inherited, no fork change) | Shared library does not build on Cygwin |
| 38548 | 2013-03-19 | - | Accuracy problem | none | not reviewed | Rounding issues in gsl-histogram with integer numbers |
| 39056 | 2013-05-23 | bug | Runtime error | none | **fixed** — Monajemi case (e4c4ac326) and the Wolpert vector now enabled (2aa89bae8) | gsl_sf_hyperg_2F1_e fails for some test cases |
| 39057 | 2013-05-23 | bug | Runtime error | none | **fixed** — gamma inverse reworked (9befdae95); report's expected value was the forward CDF | gsl_cdf_chisq_Pinv fails for some values |
| 39120 | 2013-05-29 | - | Build | none | not reviewed | Possible removal of some files |
| 39152 | 2013-06-03 | - | Accuracy problem | none | rejected — icc-only; modules pass on MSVC and gcc | make check errors with Intel icc 13.0.1 |
| 39165 | 2013-06-04 | - | Build | none | not reviewed | conditional arguments always evaluating to true |
| 39171 | 2013-06-05 | - | Performance | none | rejected — not-a-bug: GSL does not support -ffast-math | make check errors with gcc -ffast-math (or default Intel icc) |
| 39292 | 2013-06-19 | bug | Runtime error | inline | **fixed** — coulomb C = 0.5 sqrt(1+4Q) + guard vector (90d9037a5) | possible error in gsl_sf_coulomb_wave_FG_e |
| 39372 | 2013-06-30 | bug | Runtime error | none | **fixed** — gsl_hypot3 returns +Inf for any infinite argument (c1df353ae) | add check for inf/nan in gsl_hypot3 |
| 39473 | 2013-07-12 | - | Performance | clean coupling3j.patch | **fixed** — 3j symbol by edge recursion (1b2c8ff98) | more efficient algorithm for 3j,6j,9j calculations (gsl_sf_coupling_{3j,6j,9j}_e |
| 39713 | 2013-08-07 | - | Runtime error | source gsl-secant.c +3 | not reviewed | roots/secant.c "derivative value is not finite" for a good guess |
| 40092 | 2013-09-23 | - | Performance | source gsl-falsepos64.c | not reviewed | false position root finding requires too many function evals |
| 40116 | 2013-09-26 | bug | Runtime error | none | not reviewed | possible error in integration routines |
| 40176 | 2013-10-04 | bug | Runtime error | none | not reviewed | possible error in poly test suite |
| 40196 | 2013-10-07 | - | Documentation | none | **fixed** — the out-of-range key coercion (< 1 -> GAUSS15, > 6 -> GAUSS61) documented (6feb963d3) | Document gsl_integration_qag behavior on key out-of-range |
| 40755 | 2013-11-30 | bug | Accuracy problem | none | **fixed** — double cast in the Jn/Yn asymptotics test (a43fc0055) | Sporadic nan's from gsl_sf_bessel_Jn an related functions |
| 41457 | 2014-02-04 | - | Runtime error | none | not reviewed | valgrind finds errors in matrix/test.c |
| 41527 | 2014-02-09 | bug | - | none | not reviewed | Change/add multimin functions to return error codes |
| 41605 | 2014-02-15 | - | Documentation | none | **fixed** — already complete in both manuals (histogram.rst / histogram.texi); no fork change needed | gsl_histogram_pdf docs |
| 41837 | 2014-03-11 | - | Runtime error | none | **fixed** (inherited) — the three reported values are correct on the current build; the finite-sum-skip concern does not reproduce | bugs in gsl_sf_hyperg_U |
| 42042 | 2014-04-03 | bug | Runtime error | none | **fixed** — non-vanishing half-integer Jnu endpoint (e6e34279a) | nan bug in bessel_Jnu |
| 42058 | 2014-04-05 | bug | - | none | not reviewed | GSL RSS Feed does not validate |
| 42219 | 2014-04-28 | - | Runtime error | source bug_gnewton.c | not reviewed | Division by zero in "gnewton" when "f" and "fdf" differ |
| 42220 | 2014-04-28 | - | Runtime error | source bug_hybrid.c | not reviewed | Division by zero in "hybrid*" when initial guess is root |
| 42472 | 2014-05-31 | - | Runtime error | partial gsl_hh_test.c +1 | **fixed** — HH_solve/HH_svx shared tall QR (1f84bed77) | gsl_linalg_HH_solve bugs |
| 42502 | 2014-06-03 | bug | Runtime error | none | **rejected** — reporter's program omits `gsl_cdf.h`, so the function is implicitly `int`; `Pinv(0.5) = 0.0` on the built DLL and is already tested | wrong results of the function gsl_cdf_ugaussian_Pinv |
| 42830 | 2014-07-23 | bug | - | none | not reviewed | Bug in gsl_bspline_knot constructor |
| 43256 | 2014-09-19 | - | Runtime error | source sixjsymbols.c | **fixed** — stable 6j by Schulten-Gordon recurrence (59fc479e2) | gsl_sf_coupling_6j overflows |
| 43259 | 2014-09-19 | - | - | none | not reviewed | accuracy problems in specfunc |
| 43326 | 2014-09-29 | bug | Runtime error | clean bug43326.diff | **fixed** — poisson_pdf at mu=0 (4f9f4f4fc) | Bug in gsl_ran_poisson_pdf() for mu = 0.0 |
| 43496 | 2014-10-29 | bug | Runtime error | none | **fixed** — Brent parabolic-step test (c27559ec6) | Possible error in brent minimizer convergence criteria |
| 43809 | 2014-12-12 | bug | Runtime error | source gsl_hyperg.c | **fixed** — direct series for `a < 0`, large `x` (52505315d) | bug in gsl_sf_hyperg_1F1 |
| 43902 | 2014-12-29 | - | - | clean 0001-Make-the-vector-write-example-consistent-with-the-ve.patch | **fixed** — vectorw example writes 10 elements (0e82d8223) | Make the vector write example consistent with the vector read example. |
| 44612 | 2015-03-23 | bug | Runtime error | source gsl_vegas_bug_demo.c | **fixed** — VEGAS inf weight for subnormal variance (e3cab2d02) | Bug in vegas.c |
| 44865 | 2015-04-17 | - | Accuracy problem | inline | rejected — patch shortens the interval; FMA failure not reproducible here | bsimp/msbdf e5_bigt in ode-initval2/test.c is FMA-sensitive |
| 44952 | 2015-04-28 | - | Build | clean stdarg.patch | **fixed** — test/results.c uses stdarg.h (e2354de55) | test/results.c uses outdated header file name varargs.h |
| 45053 | 2015-05-07 | - | Runtime error | source gsl_bug.c +1 | rejected — bracketing patch collapses interval; 3 min tests fail | gsl_min_find_bracket is most likely incorrectly implemented |
| 45099 | 2015-05-13 | bug | - | none | not reviewed | wrong formula for BFGS update in gsl_multimin? |
| 45234 | 2015-06-02 | doc | Documentation | none | **fixed** — specfunc-mathieu.{rst,texi} already use the SF API (_e variants, int returns); no fork change needed | Mathieu function documentation hasn't been updated after switching to SF API conventions |
| 45265 | 2015-06-06 | bug | Accuracy problem | none | rejected — not reproducible on MSVC x64 (true/err ≤ 0.57 over x in [4,1000]) | gsl_sf_bessel_J0_e underestimates error for x>4 |
| 45726 | 2015-08-10 | bug | Accuracy problem | source gslbesselytest.c | **fixed** — Y family libm sin/cos (a08ef2f7f) | Incorrect results of functions bessel_y0, y1 and y2 |
| 45746 | 2015-08-13 | bug | Accuracy problem | source gsltrigtest.c | **fixed** - exact argument reduction (8a46ec7cf, 5639c380f) | Incorrect results of trigonometric functions gsl_sf_sin and gsl_sf_cos |
| 45782 | 2015-08-17 | feature | Accuracy problem | none | not reviewed | Feature request: Make derivative epsilon configurable |
| 45797 | 2015-08-19 | - | Accuracy problem | none | not reviewed | Possible problem with LAPACK Fortran routine ZGESVD |
| 45924 | 2015-09-11 | bug | Runtime error | none | **fixed** — beta inverse reworked around `t = logit(x)` with a bracketed solver (51a63cbd5) | Bug in the inverse beta function gsl_cdf_beta_Pinv, and suggested fix |
| 45925 | 2015-09-11 | - | Runtime error | none | **rejected** — not-a-bug: report confused Gamma(a,x) with P; Q matches (3cc1aaa54 docs) | Incomplete Gamma Functions flipped? |
| 46593 | 2015-12-02 | bug | Accuracy problem | none | rejected — 32-bit multifit; passes on x64 | multifit test failure in 32 bit mode |
| 46677 | 2015-12-12 | - | - | none | **deferred** — feature request: port the Wigner d-matrix (gsl_sf_wigner_drot) from contrib/wigner.c; new API, out of scope | Wigner d-matrix |
| 46678 | 2015-12-12 | bug | - | none | **fixed** — F_array spike was the #39292 Steed C; vector added (861de2adc) | Bug in gsl_sf_coulomb_wave_F_array |
| 47027 | 2016-01-31 | bug | - | none | **fixed** — WKB G' overflow now reported (8031f18bf) | gsl_sf_coulomb_wave_FG_e returns NaN but with success flag |
| 47028 | 2016-01-31 | bug | Runtime error | none | rejected — ppc64le multifit; passes on x64 | multifit testsuite failure on ppc64le |
| 47193 | 2016-02-18 | - | - | none | not reviewed | gsl_ran_poisson_pdf with mu=0 |
| 47345 | 2016-03-05 | bug | Accuracy problem | partial gsl_complex_arccosh.diff +1 | **fixed** — complex arccosh returns +0 (d6ec47d87) | arccosh(1) wrong sign |
| 47348 | 2016-03-05 | bug | - | none | not reviewed | Use of incorrect ideom floor(x+0.5) |
| 47402 | 2016-03-13 | - | - | none | not reviewed | Mathieu functions |
| 47646 | 2016-04-07 | bug | Accuracy problem | partial test_beta_small.c +2 | **fixed** — code fix upstream (05c5b5179); regression test added (d64cc4d93) | gsl_ran_beta returns NaN for small arguments |
| 48702 | 2016-08-04 | - | Runtime error | none | **fixed** — same NaN rejection as #31362 (14d595eb8) | gsl_sf_ellint_Kcomp stalls on GSL_NAN |
| 48915 | 2016-08-26 | - | Runtime error | none | rejected — AIX-only; modules pass on x64 | some test failures on AIX system for GSL 2.1.91 |
| 49465 | 2016-10-28 | - | Performance | clean 0001-initialize-newton-steffenson-solvers-with-GSL_FN_FDF.patch | **fixed** — roots Newton/Steffenson fdf init (4eac6584f) | initialize newton and steffenson solvers with GSL_FN_FDF_EVAL_F_DF |
| 49518 | 2016-11-02 | bug | Build | none | rejected — MSVC `_mktemp` crash already fixed upstream | bug in matrix/vector tests |
| 49697 | 2016-11-24 | bug | - | none | rejected — gcc -mavx passes on gcc 15.2 | gsl 2.2.1 linalg test fails with gcc (4.9.4 and later) and -mavx |
| 50343 | 2017-02-17 | - | Runtime error | inline | rejected — not-a-bug: Mathieu branch cut is arbitrary, GSL correct | Different value for mathieu_ce in Mathematica and GSL |
| 50382 | 2017-02-22 | feature | Build | none | not reviewed | Add CMAKE and NUGET support for Windows |
| 50459 | 2017-03-04 | - | - | none | **fixed** — negative-a recurrence guard for abs(a) > 2^53 (932087bcd) | Non termination of the incomplete gamma function due to floating-point rounding errors |
| 50711 | 2017-04-03 | - | - | none | **fixed** — terminating 2F1 for a non-positive-integer order (2aa89bae8) | Gauss hypergeometric function : gsl_sf_hyperg_2F1 gives up (GSL_EUNIMPL) |
| 50712 | 2017-04-03 | bug | - | none | not reviewed | Test failure for lm+accel and fdfvv |
| 50734 | 2017-04-05 | - | Performance | clean 0001-initialize-newton-steffenson-solvers-with-GSL_FN_FDF.patch | **fixed** — same fix as #49465 (4eac6584f) | initialize newton, steffenson solvers with GSL_FN_FDF_EVAL_F_DF |
| 51000 | 2017-05-11 | bug | Accuracy problem | none | **fixed** — airy_deriv at huge arguments (ac5f72d98) | Incorrect results of gsl_airy_deriv function |
| 51104 | 2017-05-24 | - | Performance | none | not reviewed | gsl_permutation_next efficiency |
| 52127 | 2017-09-27 | - | Accuracy problem | source nonsymm.c | rejected — 32 vs 64-bit eigen; x64 gives the documented answer | Difference between 32- vs. 64-bit versions of gsl_eigen_nonsymm |
| 52321 | 2017-11-01 | bug | Runtime error | clean bidiag.c.patch | rejected — patch breaks working bidiag_unpack2 | gsl_linalg_bidiag_unpack2 functioan has wrong householder transform call for V in GSL1.8 |
| 52322 | 2017-11-01 | - | Runtime error | none | rejected — 32 vs 64-bit multifit; passes on x64 | gsl_multifit_linear's output differs on 32 bit vs 64 bit linux |
| 52351 | 2017-11-06 | - | Runtime error | none | not reviewed | akima.c array indexing |
| 52359 | 2017-11-07 | - | Runtime error | source airy_divbyzero.c | **fixed** — Airy err divided by vanishing series (ac5f72d98) | Unexpected results in airy_Ai function |
| 52570 | 2017-12-01 | - | Runtime error | none | **fixed** — Airy huge negative arguments (ac5f72d98) | Inaccuracy of the Airy function due to invocation of GSL's cosine function with large input parameters |
| 52927 | 2018-01-18 | bug | Runtime error | none | rejected — not reproducible; j2 large-x tests disabled under #45730 | make check fails on Bessel j2 test |
| 53451 | 2018-03-24 | bug | - | none | **fixed** — Cauchy principal value for elliptic Pi/RJ (2673ce0d8) | gsl_sf_ellint_Pcomp( k, n, mode ) returns NaN if mode < -1 |
| 53876 | 2018-05-11 | - | Accuracy problem | none | **fixed** — missing x^(1-c) factor in the 2F1 renorm functions (46b7c412e) | gsl_sf_hyperg_2F1_renorm missing factor |
| 53903 | 2018-05-14 | bug | Runtime error | none | not reviewed | Test failure with gsl_sf_synchrotron_1_e on x86 |
| 53904 | 2018-05-14 | bug | - | none | not reviewed | Bug gsl_matrix_complex_set |
| 53905 | 2018-05-14 | bug | - | none | **fixed** — terminating 2F1 for a non-positive-integer order (2aa89bae8) | Bug in Hypergeometric function |
| 53919 | 2018-05-16 | - | Runtime error | clean v2-erf.diff +1 | **fixed** — erfc/log_erfc overflow rewritten (708791c25) | handle large values correctly in (log_)erf(c) functions |
| 54077 | 2018-06-07 | - | Runtime error | clean 0001-replace-atol-by-strtoul-in-gsl-randist.c.patch | **fixed** — gsl-randist seed via strtoul (aebe57a5c) | usage of atol in gsl-randistdoes not allow to pass big seed |
| 54919 | 2018-10-30 | bug | Build | none | rejected — icc-only; modules pass on MSVC and gcc | gsl 2.5+ test fails with icc (2016.4 and later) |
| 54925 | 2018-10-31 | - | Performance | clean 0001-Reduce-cache-misses-for-source_gemm_r.patch | rejected — out-of-scope: performance-only loop reorder | Reduce cache misses for source_gemm_r |
| 54998 | 2018-11-10 | - | Accuracy problem | none | **fixed** — same 2F1 integer-d fix (e4c4ac326, 882c8361d) | Bugs in gsl_sf_hyperg_2F1 |
| 55687 | 2019-02-10 | bug | - | none | **fixed** — NaN propagates; the `b = NaN` recursion crashed (99a73dd36) | Bad error handling in gsl_sf_hyperg_1F1_e with NaN arguments |
| 55965 | 2019-03-20 | feature | - | none | not reviewed | Implement PCG random number generator |
| 56843 | 2019-08-31 | - | Accuracy problem | none | rejected — non-x86 eigen accuracy; x64 passes | Unit Tests in linalg eigen fail on non-x86 hardware due to slight accuracy differences |
| 57173 | 2019-11-05 | feature | Accuracy problem | none | not reviewed | Feature request: zeta function for complex arguments |
| 57978 | 2020-03-09 | bug | Accuracy problem | clean sincos_pi.c.patch | **fixed** — sin_pi/cos_pi(inf) -> EDOM (7c27b358b) | Incorrect result from cosine function with inf input |
| 57979 | 2020-03-09 | bug | - | none | **fixed** — gsl_sf_hypot: +Inf for infinite args, NaN propagation (1a470222a) | Incorrect Result from hypot function with NaN input |
| 58031 | 2020-03-23 | - | Accuracy problem | none | rejected — not-a-bug: `e^x K_0(x) -> 0` is the limit; value+vector added (fc4059b12) | gsl_sf_bessel_Kn_scaled incorrectly evaluating limit |
| 58032 | 2020-03-23 | bug | Accuracy problem | none | **fixed** — 1F1 poles at x=0 (1f2607a0b) | gsl_sf_hyperg_1F1 returning incorrect result for nonpositive integers |
| 58060 | 2020-03-28 | bug | Accuracy problem | none | rejected — not-a-bug: `L^a_0 = 1`; documented + NaN rule (fc4059b12) | Dropped NaN from gsl_sf_laguerre_n |
| 58061 | 2020-03-28 | bug | Accuracy problem | none | rejected — not-a-bug: `U(0,1,x) = 1`; documented + NaN rule (fc4059b12) | Dropped NaN from gsl_sf_hyperg_U_int |
| 58062 | 2020-03-28 | bug | Accuracy problem | none | rejected — not-a-bug: `1F1(0,1,x) = 1`; reporter conceded (fc4059b12) | Dropped NaN from gsl_sf_hyperg1F1_int |
| 58063 | 2020-03-28 | - | Accuracy problem | none | rejected — not-a-bug: `hermite`/`gegenpoly` n=0 values correct; `bessel_zero_Jnu` s=0 documented (fc4059b12) | Dropped NaNs in multiple cases |
| 58064 | 2020-03-28 | - | Accuracy problem | none | **fixed** — Si/Ci limits at inf; the rest already correct (3afec212b) | Inconsistent evaluation of limits in special functions |
| 58065 | 2020-03-28 | - | - | none | **fixed** — `C_n^(0)(x) = 0` for n >= 1 (01d7e7609) | gsl_sf_gegenpoly_n returning unexpected values |
| 58066 | 2020-03-28 | bug | Accuracy problem | clean psi.c.patch | **fixed** — psi poles at all non-positive integers (d53d00109) | Digamma function returning incorrect values for (most) negative integers |
| 58067 | 2020-03-28 | - | Accuracy problem | none | rejected — not-a-bug: Ai(113) genuinely underflows; documented | Missing asymptotic behavior of the airy Ai function |
| 58068 | 2020-03-28 | - | Documentation | none | **fixed** — domain x > 0 (including x = 0 rejected) documented (42f33104b) | gsl_sf_bessel_Jnu docs lacking important domain information |
| 58069 | 2020-03-28 | doc | Documentation | none | **fixed** — Q documented for a >= 0 with Q(0,x) = 0; P deliberately left at a > 0 (it rejects a = 0) (e90d4f002) | Correction to gsl_sf_gamma_inc_Q documentation |
| 58763 | 2020-07-14 | bug | Runtime error | source demo.c +2 | rejected — not reproducible; `brent_init` already sets `c`/`fc` | gsl_root_fsolver_bren produces wrong results when run under valgrind |
| 59759 | 2020-12-23 | bug | Runtime error | none | rejected — report is about 32-bit; spmatrix passes on x64 | spmatrix test fails on x86_64 |
| 59834 | 2021-01-06 | bug | Runtime error | clean mmacc.c.patch | **fixed** — movstat accumulator alignment (4e4a88242) | Misaligned memory access error in deque.c |
| 59845 | 2021-01-08 | - | Performance | clean 0001-initialize-newton-steffenson-solvers-with-GSL_FN_FDF.patch | **fixed** — same fix as #49465 (4eac6584f) | Initialize newton, steffenson solvers with GSL_FN_FDF_EVAL_F_DF |
| 59900 | 2021-01-17 | feature | - | inline | not reviewed | Add truncated normal distribution |
| 59911 | 2021-01-20 | - | - | none | not reviewed | Problem with qagui 1D integrator |
| 59912 | 2021-01-20 | - | Documentation | none | **fixed** — gsl_permute/vector/matrix header files named in the permutation chapter (aa796e20b) | gsl_permutation header files |
| 59913 | 2021-01-20 | - | - | none | not reviewed | gsl 2.3.0 problem in gsl_integration_cquad |
| 59914 | 2021-01-20 | - | - | none | not reviewed | Native build of GSL-2.5 on windows 10 |
| 60026 | 2021-02-09 | - | - | none | not reviewed | Incorporate MIXMAX random number extension into GSL |
| 60371 | 2021-04-11 | bug | Runtime error | partial bug_interp2d_domain_error_handling.c +1 | **fixed** — interp2d domain error writes NaN (2b4e2f1e5) | Interpolation domain error handling |
| 60457 | 2021-04-26 | feature | - | clean 0001-ignore-test-files-and-doc-examples.patch +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature request: complex tridiagonal solvers |
| 60635 | 2021-05-19 | - | Accuracy problem | inline | **fixed** — const updated to CODATA 2022 (db695c8f2) | physical constants may need updating |
| 60741 | 2021-06-07 | bug | Runtime error | clean bug_60741.patch | **fixed** — lambert_W0 small arguments (fb2b4e4cd) | Inaccurate Results for Lambert W function |
| 61342 | 2021-10-16 | bug | Runtime error | inline | rejected — not-a-bug: test_c11 under -ffast-math | test_c11 test fails on ppc64 and sparc |
| 63519 | 2022-12-13 | - | Performance | none | not reviewed | gsl_root_fsolver_set "endpoints do not straddle y=0" |
| 63927 | 2023-03-14 | feature | Build | none | not reviewed | Please add a gsl-without-cblas.pc |
| 64549 | 2023-08-12 | test | Performance | clean testcases.patch | deferred — interpolation tests: leaks/dupes, low value as posted | diff patch containing new testcases for the interpolation module |
| 64613 | 2023-08-30 | - | Accuracy problem | inline | **fixed** — cdf beta_inc: 3 defects under fp-contract (1b1d94ee9) | fp-contract=fast stops convergence in beta_inc_AXPY/beta_cont_frac |
| 64777 | 2023-10-14 | bug | Runtime error | inline | rejected — not-a-bug: reporter's matrix is singular | gsl_linalg_complex_LU_decomp returns incorrect results |
| 64851 | 2023-11-03 | bug | - | clean inline gsl-config.in.patch | superseded — gsl-config exit status already done in the fork | gsl-config does not set correct status code |
| 65728 | 2024-05-12 | feature | - | none | not reviewed | Add sparse functionalities |
| 65760 | 2024-05-19 | bug | Accuracy problem | clean erfc.c.patch | **fixed** — erfc/log_erfc overflow rewritten (708791c25) | gsl_sf_log_erfc (and gsl_sf_erc) return NaN for infinite or very large finite arguments |
| 65868 | 2024-06-11 | - | - | clean 0001-bspline-Add-missing-definition-for-function.patch | **fixed** — gsl_bspline_eval_nonzero declared (ef3940132) | Missing definition for gsl_bspline_eval_nonzero |
| 65912 | 2024-06-23 | bug | - | clean 0001-Correct-GSL_SET_COMPLEX-if-native-complex-available.patch +1 | **fixed** — GSL_SET_COMPLEX signed zero / non-finite (631da98f8) | GSL_SET_COMPLEX is wrong if complex.h has no support for imaginary numbers |
| 65932 | 2024-06-30 | - | - | none | not reviewed | Complex tridiagonal solvers - patch ignored? |
| 66026 | 2024-07-26 | - | Runtime error | partial bug_66026.patch +1 | **fixed** — LU_decomp_L3 identity-initialises ipiv (d73ba5300) | gsl_linalg_LU_decomp using uninitialized memory |
| 66128 | 2024-08-27 | doc | Documentation | clean specfunc.rst.patch | **fixed** — specfunc.rst result/val typo (066f5b747) | documentation of type gsl_sf_result: confused result and val |
| 66573 | 2024-12-18 | feature | - | clean vector_complex_conjugate.patch | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: gsl_vector_complex_conjugate() |
| 66574 | 2024-12-18 | feature | - | clean 2_givens_cmplx.patch +3 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: complex singular value decomposition |
| 66575 | 2024-12-18 | feature | - | clean 5_svd_SV_solve_cmplx.patch +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: Extension to complex SVD |
| 66576 | 2024-12-18 | feature | - | clean fsolver_set_with_values.patch | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: gsl_root_fsolver_set_with_values() |
| 66695 | 2025-01-22 | feature | - | clean feagin_verner.patch | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: additional high order ODE solvers (Feagin, Verner) |
| 66742 | 2025-01-31 | - | Documentation | clean gams.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | GAMS classification |
| 66767 | 2025-02-08 | feature | Build | source binomialinv.c +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: inverse of binomial distribution |
| 66775 | 2025-02-10 | feature | Build | source poissoninv.c +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: inverse poisson distribution |
| 66800 | 2025-02-15 | feature | Accuracy problem | partial nbinomial-2.diff +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: gsl_ran_negative_binomial_pdf with p = 1 |
| 66808 | 2025-02-17 | bug | Accuracy problem | source test_airy_zeroes_derivs.c | **partial** — Airy accuracy not improved; sub-defects fixed (ac5f72d98) | Bug: Airy Ai function values inaccurate |
| 66816 | 2025-02-19 | feature | Build | source nakagami_.c | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: Nakagami random distribution |
| 66826 | 2025-02-21 | feature | Build | none | not reviewed | Feature: test cases for function hyperg_1F1() |
| 66834 | 2025-02-23 | feature | Build | clean specfunc_gamma_test.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: gamma_inc(0, 0) handling and test cases |
| 66842 | 2025-02-24 | feature | Build | clean expint_infinity.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: exponential integrals at origin |
| 66844 | 2025-02-25 | feature | Build | none | not reviewed | Feature: test cases for dilogarithm function |
| 66849 | 2025-02-26 | - | - | none | not reviewed | gsl_odeiv2_evolve_apply() may exceed final time |
| 66850 | 2025-02-26 | feature | Build | clean specfunc_hyperg_2F1.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: special cases hypergeometric2F1(a+1, b, a, x) and hypergeometric2F1(a, b+1, b, x) |
| 66862 | 2025-03-02 | bug | - | partial specfunc_test_sf.diff +1 | rejected — not-a-bug: duplicate complex sin/cos; submitter withdrew | Bug: duplicate functions in GSL produce differing results |
| 66874 | 2025-03-05 | feature | - | source doc_bst.rst | not reviewed | Feature: documentation for binary search trees (bst module) |
| 66877 | 2025-03-06 | feature | Build | none | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: test cases for hyperg_0F1() |
| 66880 | 2025-03-07 | feature | Build | clean test_hyperg_U.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: test cases for function hyperg_0F1() |
| 66886 | 2025-03-09 | feature | Build | clean specfunc_trig_.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: refactoring of specfunc/trig.c |
| 66894 | 2025-03-11 | - | Documentation | none | **fixed** — chapter introduction notes that randist does not validate distribution parameters (dbe316d49) | Domain value checking for random number distributions |
| 66922 | 2025-03-17 | feature | Build | clean specfunc_trig.diff | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: special case for special functin logsin() |
| 66949 | 2025-03-25 | feature | Build | source erlang.c | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: Erlang cumulative distribution |
| 66993 | 2025-04-05 | bug | Accuracy problem | none | not reviewed | Bug: pow_int issues |
| 67058 | 2025-04-28 | - | Accuracy problem | none | **fixed** — empty data sets raise `GSL_EBADLEN` (NaN with the handler off) (68fc3c752) | gsl_stats_mean and gsl_stats_sd result to 0.0 when array with zero length |
| 67301 | 2025-07-10 | bug | - | none | not reviewed | Bug: Test for existence of uniform random variate in histogram |
| 67359 | 2025-07-27 | feature | - | clean hyperg_2F0.c.patch +1 | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: Add special case for gsl_sf_hyperg_2F0 |
| 67445 | 2025-08-20 | - | Build | clean 0001-linalg-increase-cholesky_invert-Hilbert-test-toleran.patch | **fixed** — cholesky_invert Hilbert tolerance (f84a57a0f) | linalg test failures under gcc 14.2.1 |
| 67446 | 2025-08-20 | - | Build | none | rejected — gcc 14.2.1; passes on gcc 15.2 x64 | multilarge_nlinear test failures under gcc 14.2.1 |
| 67447 | 2025-08-20 | - | Build | none | rejected — gcc 14.2.1; passes on gcc 15.2 x64 | spmatrix test failures under gcc 14.2.1 |
| 67494 | 2025-09-09 | - | Accuracy problem | none | not reviewed | Consistency in pow_int usage |
| 67621 | 2025-10-23 | - | Documentation | clean rst_deprecated.diff | **fixed** — deprecation markers from headers (5df6c0079) | Mark functions as deprecated if they are also deprecated in the code |
| 67689 | 2025-11-10 | doc | Documentation | clean specfunc-psi.rst.patch | **fixed** — gsl_sf_complex_psi_e documented (345883172) | Incomplete documentation of digamma functions in GSL specfunc |
| 67705 | 2025-11-15 | bug | Build | none | rejected — not reproducible; linalg passes with a fixed seed | ttest failure in linalg/QR_solve_r random |
| 67728 | 2025-11-23 | bug | Accuracy problem | source psi_dropin.c | **fixed** — polygamma at negative arguments (e9f69933a) | gsl_sf_psi_n_e yields domain error |
| 67774 | 2025-12-05 | feature | Accuracy problem | none | not reviewed | Feature: arctan integral is also defined for negative inputs |
| 68068 | 2026-02-19 | bug | Accuracy problem | none | not reviewed | Bug: incorrect straddling of area of convergence in quad_golden |
| 68073 | 2026-02-20 | bug | Documentation | none | **fixed** — gsl_stats_select comment corrected from "k-th largest" to "k-th smallest" (521bed61d) | Bug: incorrect inline code comment on BASE FUNCTION(gsl_stats,select) |
| 68098 | 2026-02-27 | feature | - | none | not reviewed | Feature: division by zero when data is perfectly correlated |
| 68283 | 2026-04-26 | - | - | none | not reviewed | Correction to gsl_rstat_skew and gsl_rstat_kurtosis |
| 68312 | 2026-05-07 | bug | Accuracy problem | none | **fixed** — same recurrence as #43256 (59fc479e2) | Wigner symbols inaccurate for large j |
| 68367 | 2026-05-19 | feature | - | source invelljac.c | not reviewed (feature-shaped) — skipped in first pass as feature-shaped | Feature: inverse Jacobi elliptic integrals |
| 68379 | 2026-05-21 | - | - | none | not reviewed | Histogram: expand scope of internal variables |
| 68398 | 2026-05-26 | feature | - | none | not reviewed | Feature: use GSL native gsl_ldexp in eta fuction for integer argument |
| 68415 | 2026-06-02 | bug | - | none | not reviewed | Bug: inconsistency in variance error handling in statictics module |
| 68479 | 2026-06-25 | bug | Accuracy problem | none | not reviewed | Bug: gsl_ran_binomial is not accurate |
| 68495 | 2026-07-03 | bug | Build | inline | **fixed** — pow_int INT_MIN overflow (41b1e2c00) | Undefined behavior in gsl_pow_int: signed int overflow |
| 68518 | 2026-07-13 | - | Build | none | not reviewed | Remove stale construct in configure.ac |
| 68549 | 2026-07-23 | - | - | none | not reviewed | Breaking change: match argument list among distribution functions |
| 68592 | 2026-08-03 | doc | Documentation | none | **fixed** — new BST chapter added to both manuals and to the index/menu, marked as an unreviewed AI draft (20eda0316) | Documentation: missing Binary search tree documentation from index |
| 68611 | 2026-08-13 | feature | - | none | not reviewed | Feature: guard against degenerate histogram input |
| 68625 | 2026-08-21 | bug | Accuracy problem | none | **fixed** — hermite_func_der at n=0,1 (5d36e1999) | Bug: incorrect derivative of the Hermite function of order 0 or 1 |
| 68663 | 2026-08-31 | feature | - | none | not reviewed | Feature: remove GAMM_INC macro in favor of direct call |
| 68704 | 2026-09-18 | feature | - | none | not reviewed | Feature: introduce variance of a histogram gsl_histogram_variance |
