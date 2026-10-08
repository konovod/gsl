# Savannah open-issue triage

Generated 2026-10-02 from the `scripts/savannah_bugs.py` sweep (`inventory.tsv`).

Scratch index, not part of the fork's record of changes. Verdicts are curated from
`FORKNEWS`, `SAVANNAH_REVIEW.md` and the git history; every other row is un-triaged.

**219 open items: 219 triaged, 0 never reviewed — 191 closed or rejected, 28 parked in the backlog below.**

| status | count |
|---|---:|
| fixed | 150 |
| partial | 0 |
| rejected | 40 |
| backlog | 28 |
| superseded | 1 |
| not reviewed | 0 |

Legend for the *Patch* column: the sweep's `git apply` verdict (`clean` / `partial` / `dirty`), `source` for a whole
replacement file, `inline` for a patch pasted into the bug text, `none` for no patch. The **backlog** verdict marks a
report whose reported defect is closed or declined but which still carries work, listed under *Future work backlog*
below; it is not counted as rejected.

| # | Date | Kind | Cat | Patch | Fork status | Notes |
|---|---|---|---|---|---|---|
| 21828 | 2007-12-18 | - | Performance | none | **fixed** — `gsl_linalg_householder_hm` applied the reflector column-by-column on row-major storage, and the accessors were out-of-line calls without `HAVE_INLINE`; rewritten to walk rows in storage order, bit-identical (cd9e31326) | suboptimal performance of gsl_fdfsolver_lmsder |
| 21831 | 2007-12-18 | - | Accuracy problem | source levy.c | **rejected** — not a bug: the CMS transform is exact for alpha < 1 (characteristic-function check); `beta = 0` delegates to `gsl_ran_levy` | Levý random number generator for alpha < 1 |
| 21833 | 2007-12-18 | - | Performance | none | **fixed** — the selection loop scans backward from the suffix (c382ed86e); the alleged quadratic behaviour is amortised O(1) per permutation and the 8! walk in `permutation/test.c` covers it | suboptimal performance of gsl permutation? |
| 21835 | 2007-12-18 | - | Accuracy problem | clean test_hyperg.diff +1 | **backlog** — non-positive-integer termination for `x >= 1` fixed (2aa89bae8); the `c = a+b` near `x = 1` case still returns `GSL_EMAXITER` (residual defect below) | gsl_sf_hyperg_2F1 problematic arguments |
| 21836 | 2007-12-18 | - | Accuracy problem | none | **fixed** — Q returns exact complement of P in the series window (6f23a4cb5) | gamma_inc_P and gamma_inc_Q only satisfy P+Q=1 within errors |
| 21837 | 2007-12-18 | - | Runtime error | none | **backlog** — docs corrected (9418afa5f); a zero-diagonal permutation path is new algorithm work | gsl_linalg_solve_symm_tridiag requires positive definite matrix |
| 24252 | 2008-09-12 | feature | - | source gamma_tail_jpl_080908.c | **backlog** — the posted truncated-gamma sampler is wrong (the `a<1` branch uses `a` instead of the tail); a correct gamma-tail generator is wanted work; see FORKNEWS | suggestion: add gamma tail distribution |
| 24871 | 2008-11-18 | feature | - | none | **fixed** — complex exponential integrals E1(z), En(z), their exp(z)-scaled variants and Ei(z), by a port of Amos Algorithm 683 (c01e6b922, 4416fe729, ffad9c422) | suggestion, add support for E_n |
| 25320 | 2009-01-14 | - | Accuracy problem | none | **backlog** — no Fresnel in the tree; importing it is new API/algorithm work | Import fresnel, bugs on GSL Extension Fresnel |
| 28267 | 2009-12-11 | - | Accuracy problem | source hyperg1F1.c | **fixed** — same defect as #43809 (52505315d); the transition region `x ~ a^2` is handled by a double-double series (135c6eb27) | poor convergence region for gsl_sf_hyperg_1F1 |
| 29834 | 2010-05-09 | - | Runtime error | source error_cblas_v2.h | **backlog** — the report's cblas checking macros already live in `cblas/error_cblas*.h` (inherited); unifying the CBLAS/GSL argument-checking macros is wanted refactoring | insufficient argument checking in blas wrapper |
| 30324 | 2010-07-02 | - | Accuracy problem | none | **fixed** — 2F1 continued to x < -1 by the Pfaff transformations, real and conjugate (b99ec4e7c, cfe9311cd) | improve range of 2F1 |
| 30510 | 2010-07-21 | - | Runtime error | none | **fixed** — U(a,b,x) for x < 0, integer b (00859f816) and non-integer b (77388cb67), on the real branch | problems with hyperg_U(a,b,x) for x<0 |
| 30540 | 2010-07-24 | feature | Accuracy problem | partial bug-ode2.c +2 | **fixed** — v1 rk2imp/rk4imp iterate to convergence; a non-converged step is rejected through yerr (fa1622e11) | please, add convergence checks in rk4imp/rk2imp |
| 30583 | 2010-07-28 | doc | Documentation | none | **fixed** — Legendre/Carlson relations and the negative-parameter (imaginary-modulus) transformation documented (b0eec8bc6) | improve documentation for Elliptic functions |
| 30885 | 2010-08-27 | - | Runtime error | none | **fixed** — Coulomb F recurrence rescaled, no overflow (25841970e) | nans from gsl_sf_coulomb_wave_FG_e(1.2693881947287221e-07, 0.0, lam_F=37, lam_G=36) |
| 30947 | 2010-09-02 | - | - | clean 0001-Fixed-step-size-control-object.patch | **fixed** — `gsl_odeiv_control_fixed_new` for the legacy v1 interface (379750b93) | Please, include fixed step size control object for ode suite |
| 31109 | 2010-09-23 | - | Performance | none | **fixed** — bsimp derives its order from the driver control's error level; loose tolerances drop from order 12 to 8 (1fbb40ca3) | ode-initval/bsimp is always high order |
| 31362 | 2010-10-18 | bug | Runtime error | none | **fixed** — NaN rejected as GSL_EDOM in the complete elliptic integrals (14d595eb8) | The Complete Elliptic Integrals (gsl_sf_ellint_Ecomp and _Kcomp) Loop Forever with NaN Argument |
| 31426 | 2010-10-23 | - | Runtime error | none | **fixed** — rescaled/bounded symmetric QR iteration (7d586d89b) | infinite loop in gsl_eigen_symm |
| 32257 | 2011-01-26 | - | - | none | **fixed** — Gauss-Lobatto, Gauss-Radau, Clenshaw-Curtis and Fejer fixed quadrature rules (7cdd25747, 4ad7f6f36) | RFE: Import integration routines from quadrule |
| 32306 | 2011-01-31 | bug | Accuracy problem | source hyp.c | **fixed** — integer-d 2F1 series and error estimate for `x < 0.995` (e4c4ac326, 882c8361d) and the integer-d reflection at `x >= 0.995` (58936d413) | sign error in gsl_sf_hyperg_2F1 |
| 32776 | 2011-03-14 | feature | - | source quadratic.c | **fixed** — multimin quadratic minimiser (53cd0d098) | RFE: Add brute-force quadratic numerical multidimensional minimizer |
| 34361 | 2011-09-22 | - | Runtime error | none | **backlog** — the `init_augment` monotonicity guard rejects non-monotone breakpoints (74cace84f; test/docs c8543330d); constrained least squares for `gsl_bspline_knots_greville` remains | gsl_bspline_knots_greville needs inequality constrained linear least squares |
| 35032 | 2011-12-11 | doc | Documentation | none | **fixed** — gsl_test support functions documented in the usage chapter, marked as an unreviewed AI draft (8ef88af6e) | gsl_test.h lacks documentation in the reference manual |
| 36152 | 2012-04-11 | bug | Accuracy problem | source testbessel.c | **backlog** — asymptotics fixed via Y-family libm sin/cos (a08ef2f7f) and exact sin/cos reduction (8a46ec7cf); a `gsl_sf_sinc_e` large-`x` test vector remains | Incorrect asymptotics of spherical Bessel functions |
| 36197 | 2012-04-15 | - | Build | dirty 36197b.diff +1 | **fixed** — leading `__` dropped from every include guard (e40dc7ceb) | reserved identifier violation |
| 36578 | 2012-06-02 | - | Documentation | none | **fixed** — doc now says the error handler is invoked with GSL_EUNSUP and that the void function returns no code (060c2472c) | inconsistency in gsl_ieee_env_setup doc/api |
| 37209 | 2012-08-28 | bug | Accuracy problem | none | **fixed** — upstream 441bc40ff (inherited); NaN gone | gsl_sf_bessel_jl_e returns NaN for large inputs without setting error code |
| 37408 | 2012-09-20 | - | - | clean au.patch | **fixed** — const AU/parsec updated to IAU 2012 (3377c7fbe) | The astronomical unit (AU) has been re-defined |
| 37894 | 2012-12-10 | bug | Build | clean gsl-autotools.diff | **fixed** — upstream 1d002ee93/ae19e3e8b (inherited, no fork change) | Shared library does not build on Cygwin |
| 38548 | 2013-03-19 | - | Accuracy problem | none | **fixed** — exact integer bin limits (f9c977f52) | Rounding issues in gsl-histogram with integer numbers |
| 39056 | 2013-05-23 | bug | Runtime error | none | **fixed** — Monajemi case (e4c4ac326) and the Wolpert vector now enabled (2aa89bae8) | gsl_sf_hyperg_2F1_e fails for some test cases |
| 39057 | 2013-05-23 | bug | Runtime error | none | **fixed** — gamma inverse reworked (9befdae95); report's expected value was the forward CDF | gsl_cdf_chisq_Pinv fails for some values |
| 39120 | 2013-05-29 | - | Build | none | **fixed** — five unreferenced files removed; modnewton1.c/matrix.c kept (14920fa04) | Possible removal of some files |
| 39152 | 2013-06-03 | - | Accuracy problem | none | **backlog** — icc-only; modules pass on MSVC and gcc; re-check parked pending Intel oneAPI | make check errors with Intel icc 13.0.1 |
| 39165 | 2013-06-04 | - | Build | none | rejected — both sites already corrected (#58065 rewrite / bzr 4829); the remaining `hmax` guard is not always true | conditional arguments always evaluating to true |
| 39171 | 2013-06-05 | - | Performance | none | rejected — not-a-bug: GSL does not support -ffast-math | make check errors with gcc -ffast-math (or default Intel icc) |
| 39292 | 2013-06-19 | bug | Runtime error | inline | **fixed** — coulomb C = 0.5 sqrt(1+4Q) + guard vector (90d9037a5) | possible error in gsl_sf_coulomb_wave_FG_e |
| 39372 | 2013-06-30 | bug | Runtime error | none | **fixed** — gsl_hypot3 returns +Inf for any infinite argument (c1df353ae) | add check for inf/nan in gsl_hypot3 |
| 39473 | 2013-07-12 | - | Performance | clean coupling3j.patch | **fixed** — 3j symbol by edge recursion (1b2c8ff98) | more efficient algorithm for 3j,6j,9j calculations (gsl_sf_coupling_{3j,6j,9j}_e |
| 39713 | 2013-08-07 | - | Runtime error | source gsl-secant.c +3 | **fixed** — inherited upstream d50dc70e5/eaaae349d; test vector added by ce90a6b91; posted patches are pre-2013 | roots/secant.c "derivative value is not finite" for a good guess |
| 40092 | 2013-09-23 | - | Performance | source gsl-falsepos64.c | **fixed** — skip the redundant linear-interpolation evaluation when it lands on an endpoint (0697ba772); the report's 64-bit function drops from 98 to 52 evaluations | false position root finding requires too many function evals |
| 40116 | 2013-09-26 | bug | Runtime error | none | rejected - faithful QUADPACK port; `large_interval`/`increase_nrmax` map to `go to 90`, no failing case | possible error in integration routines |
| 40176 | 2013-10-04 | bug | Runtime error | none | rejected — already fixed upstream (9cc12d037, 0466df866); duplicate of #39055 | possible error in poly test suite |
| 40196 | 2013-10-07 | - | Documentation | none | **fixed** — the out-of-range key coercion (< 1 -> GAUSS15, > 6 -> GAUSS61) documented (6feb963d3) | Document gsl_integration_qag behavior on key out-of-range |
| 40755 | 2013-11-30 | bug | Accuracy problem | none | **fixed** — double cast in the Jn/Yn asymptotics test (a43fc0055) | Sporadic nan's from gsl_sf_bessel_Jn an related functions |
| 41457 | 2014-02-04 | - | Runtime error | none | rejected — already fixed upstream by 6ed874986 (memset covers the padding bytes) | valgrind finds errors in matrix/test.c |
| 41527 | 2014-02-09 | bug | - | none | **fixed** — the simplex minimizers propagate GSL_EBADFUNC instead of masking a non-finite objective as GSL_EFAILED (f007aa42c) | Change/add multimin functions to return error codes |
| 41605 | 2014-02-15 | - | Documentation | none | **fixed** — already complete in both manuals (histogram.rst / histogram.texi); no fork change needed | gsl_histogram_pdf docs |
| 41837 | 2014-03-11 | - | Runtime error | none | **fixed** (inherited) — the three reported values are correct on the current build; the finite-sum-skip concern does not reproduce | bugs in gsl_sf_hyperg_U |
| 42042 | 2014-04-03 | bug | Runtime error | none | **fixed** — non-vanishing half-integer Jnu endpoint (e6e34279a) | nan bug in bessel_Jnu |
| 42058 | 2014-04-05 | bug | - | none | rejected — not a GSL defect: the feed is generated by Savannah (`news/atom.php`), not the GSL tree | GSL RSS Feed does not validate |
| 42219 | 2014-04-28 | - | Runtime error | source bug_gnewton.c | **fixed** — guard the step reduction with phi0 > 0 (8640b846e); trap-enabled regression test | Division by zero in "gnewton" when "f" and "fdf" differ |
| 42220 | 2014-04-28 | - | Runtime error | source bug_hybrid.c | **fixed** — hybrid/hybridj return at fnorm = 0 and compute_wv guards pnorm = 0 (8640b846e); trap-enabled regression test | Division by zero in "hybrid*" when initial guess is root |
| 42472 | 2014-05-31 | - | Runtime error | partial gsl_hh_test.c +1 | **fixed** — HH_solve/HH_svx shared tall QR (1f84bed77) | gsl_linalg_HH_solve bugs |
| 42502 | 2014-06-03 | bug | Runtime error | none | **rejected** — reporter's program omits `gsl_cdf.h`, so the function is implicitly `int`; `Pinv(0.5) = 0.0` on the built DLL and is already tested | wrong results of the function gsl_cdf_ugaussian_Pinv |
| 42830 | 2014-07-23 | bug | - | none | **fixed** - non-monotone breakpoints rejected in init_augment (74cace84f) | Bug in gsl_bspline_knot constructor |
| 43256 | 2014-09-19 | - | Runtime error | source sixjsymbols.c | **fixed** — stable 6j by Schulten-Gordon recurrence (59fc479e2) | gsl_sf_coupling_6j overflows |
| 43259 | 2014-09-19 | - | - | none | **backlog** — `hyperg_0F1` (5fafddc63), `exprel_2` (c2ba2e271) and `gamma_inc_Q` (18652f8c1) fixed; clausen/zeta/eta already fixed; the psi_1/pochrel/ellint_P "degenerate" inputs need a fresh verdict (residual defect below) | accuracy problems in specfunc |
| 43326 | 2014-09-29 | bug | Runtime error | clean bug43326.diff | **fixed** — poisson_pdf at mu=0 (4f9f4f4fc) | Bug in gsl_ran_poisson_pdf() for mu = 0.0 |
| 43496 | 2014-10-29 | bug | Runtime error | none | **fixed** — Brent parabolic-step test (c27559ec6) | Possible error in brent minimizer convergence criteria |
| 43809 | 2014-12-12 | bug | Runtime error | source gsl_hyperg.c | **fixed** — direct series for `a < 0`, large `x` (52505315d) | bug in gsl_sf_hyperg_1F1 |
| 43902 | 2014-12-29 | - | - | clean 0001-Make-the-vector-write-example-consistent-with-the-ve.patch | **fixed** — vectorw example writes 10 elements (0e82d8223) | Make the vector write example consistent with the vector read example. |
| 44612 | 2015-03-23 | bug | Runtime error | source gsl_vegas_bug_demo.c | **fixed** — VEGAS inf weight for subnormal variance (e3cab2d02) | Bug in vegas.c |
| 44865 | 2015-04-17 | - | Accuracy problem | inline | **backlog** — FMA failure not reproducible here and the patch shortens the interval; re-check parked pending `-ffp-contract=fast` testing | bsimp/msbdf e5_bigt in ode-initval2/test.c is FMA-sensitive |
| 44952 | 2015-04-28 | - | Build | clean stdarg.patch | **fixed** — test/results.c uses stdarg.h (e2354de55) | test/results.c uses outdated header file name varargs.h |
| 45053 | 2015-05-07 | - | Runtime error | source gsl_bug.c +1 | rejected — bracketing patch collapses interval; 3 min tests fail | gsl_min_find_bracket is most likely incorrectly implemented |
| 45099 | 2015-05-13 | bug | - | none | rejected - not-a-bug: the BFGS update is the standard dense update, H implicit in p | wrong formula for BFGS update in gsl_multimin? |
| 45234 | 2015-06-02 | doc | Documentation | none | **fixed** — specfunc-mathieu.{rst,texi} already use the SF API (_e variants, int returns); no fork change needed | Mathieu function documentation hasn't been updated after switching to SF API conventions |
| 45265 | 2015-06-06 | bug | Accuracy problem | none | rejected — not reproducible on MSVC x64 (true/err ≤ 0.57 over x in [4,1000]) | gsl_sf_bessel_J0_e underestimates error for x>4 |
| 45726 | 2015-08-10 | bug | Accuracy problem | source gslbesselytest.c | **fixed** — Y family libm sin/cos (a08ef2f7f) | Incorrect results of functions bessel_y0, y1 and y2 |
| 45746 | 2015-08-13 | bug | Accuracy problem | source gsltrigtest.c | **fixed** - exact argument reduction (8a46ec7cf, 5639c380f) | Incorrect results of trigonometric functions gsl_sf_sin and gsl_sf_cos |
| 45782 | 2015-08-17 | feature | Accuracy problem | none | **fixed** — configurable finite-difference Jacobian step on `gsl_multiroot_fsolver`, default unchanged (8abbe2c30, 296f3b58a) | Feature request: Make derivative epsilon configurable |
| 45797 | 2015-08-19 | - | Accuracy problem | none | rejected — external LAPACK/distribution report: the program links LAPACK `ZGESVD`; GSL has no complex SVD and no GSL code path exists | Possible problem with LAPACK Fortran routine ZGESVD |
| 45924 | 2015-09-11 | bug | Runtime error | none | **fixed** — beta inverse reworked around `t = logit(x)` with a bracketed solver (51a63cbd5) | Bug in the inverse beta function gsl_cdf_beta_Pinv, and suggested fix |
| 45925 | 2015-09-11 | - | Runtime error | none | **rejected** — not-a-bug: report confused Gamma(a,x) with P; Q matches (3cc1aaa54 docs) | Incomplete Gamma Functions flipped? |
| 46593 | 2015-12-02 | bug | Accuracy problem | none | rejected — 32-bit multifit; passes on x64 | multifit test failure in 32 bit mode |
| 46677 | 2015-12-12 | - | - | none | **fixed** — ported the Wigner d-matrix (`gsl_sf_wigner_drot`) into `specfunc/coupling.c` (17e39bcbd) | Wigner d-matrix |
| 46678 | 2015-12-12 | bug | - | none | **fixed** — F_array spike was the #39292 Steed C; vector added (861de2adc) | Bug in gsl_sf_coulomb_wave_F_array |
| 47027 | 2016-01-31 | bug | - | none | **fixed** — WKB G' overflow now reported (8031f18bf) | gsl_sf_coulomb_wave_FG_e returns NaN but with success flag |
| 47028 | 2016-01-31 | bug | Runtime error | none | rejected — ppc64le multifit; passes on x64 | multifit testsuite failure on ppc64le |
| 47193 | 2016-02-18 | - | - | none | **fixed** — duplicate of #43326 (4f9f4f4fc); regression test 5346bba43 | gsl_ran_poisson_pdf with mu=0 |
| 47345 | 2016-03-05 | bug | Accuracy problem | partial gsl_complex_arccosh.diff +1 | **fixed** — complex arccosh returns +0 (d6ec47d87) | arccosh(1) wrong sign |
| 47348 | 2016-03-05 | bug | - | none | **fixed** - floor(x+0.5) -> rint; hyperg_1F1 integer test also missing fabs (b681165e1); coulomb keeps deliberate round-half-up | Use of incorrect ideom floor(x+0.5) |
| 47402 | 2016-03-13 | - | - | none | **backlog** — Mathieu return types + undocumented Fourier coefficients (a31710cf6) and angular ce/se derivatives (47727ca0e) fixed; item 3 (coefficient caching) deferred, design in SAVANNAH_REVIEW.md | Mathieu functions |
| 47646 | 2016-04-07 | bug | Accuracy problem | partial test_beta_small.c +2 | **fixed** — code fix upstream (05c5b5179); regression test added (d64cc4d93) | gsl_ran_beta returns NaN for small arguments |
| 48702 | 2016-08-04 | - | Runtime error | none | **fixed** — same NaN rejection as #31362 (14d595eb8) | gsl_sf_ellint_Kcomp stalls on GSL_NAN |
| 48915 | 2016-08-26 | - | Runtime error | none | rejected — AIX-only; modules pass on x64 | some test failures on AIX system for GSL 2.1.91 |
| 49465 | 2016-10-28 | - | Performance | clean 0001-initialize-newton-steffenson-solvers-with-GSL_FN_FDF.patch | **fixed** — roots Newton/Steffenson fdf init (4eac6584f) | initialize newton and steffenson solvers with GSL_FN_FDF_EVAL_F_DF |
| 49518 | 2016-11-02 | bug | Build | none | rejected — MSVC `_mktemp` crash already fixed upstream | bug in matrix/vector tests |
| 49697 | 2016-11-24 | bug | - | none | **backlog** — gcc `-mavx` passes on gcc 15.2; re-check parked pending an explicit `-mavx` run | gsl 2.2.1 linalg test fails with gcc (4.9.4 and later) and -mavx |
| 50343 | 2017-02-17 | - | Runtime error | inline | rejected — not-a-bug: Mathieu branch cut is arbitrary, GSL correct | Different value for mathieu_ce in Mathematica and GSL |
| 50382 | 2017-02-22 | feature | Build | none | **backlog** — CMake build already in the fork; NuGet packaging is new distribution surface, wanted but out of scope of the report | Add CMAKE and NUGET support for Windows |
| 50459 | 2017-03-04 | - | - | none | **fixed** — negative-a recurrence guard for abs(a) > 2^53 (932087bcd) | Non termination of the incomplete gamma function due to floating-point rounding errors |
| 50711 | 2017-04-03 | - | - | none | **fixed** — terminating 2F1 for a non-positive-integer order (2aa89bae8) | Gauss hypergeometric function : gsl_sf_hyperg_2F1 gives up (GSL_EUNIMPL) |
| 50712 | 2017-04-03 | bug | - | none | **backlog** — reproduces on MSVC x64 (lmaccel + finite-difference fvv on box3d); the test is disabled with `#if 0` and no small fix has been found (residual defect below) | Test failure for lm+accel and fdfvv |
| 50734 | 2017-04-05 | - | Performance | clean 0001-initialize-newton-steffenson-solvers-with-GSL_FN_FDF.patch | **fixed** — same fix as #49465 (4eac6584f) | initialize newton, steffenson solvers with GSL_FN_FDF_EVAL_F_DF |
| 51000 | 2017-05-11 | bug | Accuracy problem | none | **fixed** — airy_deriv at huge arguments (ac5f72d98) | Incorrect results of gsl_airy_deriv function |
| 51104 | 2017-05-24 | - | Performance | none | **fixed** — the Algorithm L selection loop scans backward from the suffix (c382ed86e) | gsl_permutation_next efficiency |
| 52127 | 2017-09-27 | - | Accuracy problem | source nonsymm.c | rejected — 32 vs 64-bit eigen; x64 gives the documented answer | Difference between 32- vs. 64-bit versions of gsl_eigen_nonsymm |
| 52321 | 2017-11-01 | bug | Runtime error | clean bidiag.c.patch | rejected — patch breaks working bidiag_unpack2 | gsl_linalg_bidiag_unpack2 functioan has wrong householder transform call for V in GSL1.8 |
| 52322 | 2017-11-01 | - | Runtime error | none | rejected — 32 vs 64-bit multifit; passes on x64 | gsl_multifit_linear's output differs on 32 bit vs 64 bit linux |
| 52351 | 2017-11-06 | - | Runtime error | none | **fixed** — average the Akima slopes when a weight vanishes (0efb5f87f) | akima.c array indexing |
| 52359 | 2017-11-07 | - | Runtime error | source airy_divbyzero.c | **fixed** — Airy err divided by vanishing series (ac5f72d98) | Unexpected results in airy_Ai function |
| 52570 | 2017-12-01 | - | Runtime error | none | **fixed** — Airy huge negative arguments (ac5f72d98) | Inaccuracy of the Airy function due to invocation of GSL's cosine function with large input parameters |
| 52927 | 2018-01-18 | bug | Runtime error | none | rejected — not reproducible; j2 large-x tests disabled under #45730 | make check fails on Bessel j2 test |
| 53451 | 2018-03-24 | bug | - | none | **fixed** — Cauchy principal value for elliptic Pi/RJ (2673ce0d8) | gsl_sf_ellint_Pcomp( k, n, mode ) returns NaN if mode < -1 |
| 53876 | 2018-05-11 | - | Accuracy problem | none | **fixed** — missing x^(1-c) factor in the 2F1 renorm functions (46b7c412e) | gsl_sf_hyperg_2F1_renorm missing factor |
| 53903 | 2018-05-14 | bug | Runtime error | none | rejected — 32-bit x87 musl `pow` returns 0; the reported value/error are exactly the terms that survive | Test failure with gsl_sf_synchrotron_1_e on x86 |
| 53904 | 2018-05-14 | bug | - | none | rejected — not a bug: uninitialized matrix; direct-call and variable forms of `gsl_complex_rect` are equivalent | Bug gsl_matrix_complex_set |
| 53905 | 2018-05-14 | bug | - | none | **fixed** — terminating 2F1 for a non-positive-integer order (2aa89bae8) | Bug in Hypergeometric function |
| 53919 | 2018-05-16 | - | Runtime error | clean v2-erf.diff +1 | **fixed** — erfc/log_erfc overflow rewritten (708791c25); the related hazard error-bar defect also fixed (9c43a759f) | handle large values correctly in (log_)erf(c) functions |
| 54077 | 2018-06-07 | - | Runtime error | clean 0001-replace-atol-by-strtoul-in-gsl-randist.c.patch | **fixed** — gsl-randist seed via strtoul (aebe57a5c) | usage of atol in gsl-randistdoes not allow to pass big seed |
| 54919 | 2018-10-30 | bug | Build | none | **backlog** — icc-only; modules pass on MSVC and gcc; re-check parked pending Intel oneAPI | gsl 2.5+ test fails with icc (2016.4 and later) |
| 54925 | 2018-10-31 | - | Performance | clean 0001-Reduce-cache-misses-for-source_gemm_r.patch | **fixed** — loops reordered so a row of C stays resident (39cde03e2); numerically neutral | Reduce cache misses for source_gemm_r |
| 54998 | 2018-11-10 | - | Accuracy problem | none | **fixed** — same 2F1 integer-d fix (e4c4ac326, 882c8361d) | Bugs in gsl_sf_hyperg_2F1 |
| 55687 | 2019-02-10 | bug | - | none | **fixed** — NaN propagates; the `b = NaN` recursion crashed (99a73dd36) | Bad error handling in gsl_sf_hyperg_1F1_e with NaN arguments |
| 55965 | 2019-03-20 | feature | - | none | **fixed** — PCG32 added as `gsl_rng_pcg32` (3fd2c7475) | Implement PCG random number generator |
| 56843 | 2019-08-31 | - | Accuracy problem | none | **backlog** — non-x86 eigen accuracy; x64 passes; reproduce on the arm64 macOS CI target | Unit Tests in linalg eigen fail on non-x86 hardware due to slight accuracy differences |
| 57173 | 2019-11-05 | feature | Accuracy problem | none | **fixed** — complex `gsl_sf_complex_zeta_e`, `_hzeta_e` and `_eta_e` via Euler-Maclaurin and reflection (2708ca0fa, a605f0c44) | Feature request: zeta function for complex arguments |
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
| 59900 | 2021-01-17 | feature | - | inline | **fixed** — truncated Gaussian randist + cdf functions; generator self-contained in randist to avoid the cdf cycle (76f3b1db8) | Add truncated normal distribution |
| 59911 | 2021-01-20 | - | - | none | rejected — caller's cosh(x) overflows; library returns GSL_EMAXITER | Problem with qagui 1D integrator |
| 59912 | 2021-01-20 | - | Documentation | none | **fixed** — gsl_permute/vector/matrix header files named in the permutation chapter (aa796e20b) | gsl_permutation header files |
| 59913 | 2021-01-20 | - | - | none | **fixed** — cquad error estimate guarded against nc == 0 (f20496de6) | gsl 2.3.0 problem in gsl_integration_cquad |
| 59914 | 2021-01-20 | - | - | none | **backlog** — no defect: a success story offering a build write-up; a documented native Windows build write-up is wanted | Native build of GSL-2.5 on windows 10 |
| 60026 | 2021-02-09 | feature | - | none | **fixed** — MIXMAX with N=17 over 2^61-1 added as `gsl_rng_mixmax17`; CLHEP's `seed_spbox` seeding needs no 128-bit multiply and get returns the low 32 bits (260557ee9) | Incorporate MIXMAX random number extension into GSL |
| 60371 | 2021-04-11 | bug | Runtime error | partial bug_interp2d_domain_error_handling.c +1 | **fixed** — interp2d domain error writes NaN (2b4e2f1e5) | Interpolation domain error handling |
| 60457 | 2021-04-26 | feature | - | clean 0001-ignore-test-files-and-doc-examples.patch +1 | **fixed** — four complex tridiagonal solvers in `linalg/tridiagcomplex.c`, verified against numpy and the residual (dd5e0a0f2, ee6eee516) | Feature request: complex tridiagonal solvers |
| 60635 | 2021-05-19 | - | Accuracy problem | inline | **fixed** — const updated to CODATA 2022 (db695c8f2) | physical constants may need updating |
| 60741 | 2021-06-07 | bug | Runtime error | clean bug_60741.patch | **fixed** — lambert_W0 small arguments (fb2b4e4cd) | Inaccurate Results for Lambert W function |
| 61342 | 2021-10-16 | bug | Runtime error | inline | rejected — not-a-bug: test_c11 under -ffast-math | test_c11 test fails on ppc64 and sparc |
| 63519 | 2022-12-13 | - | Performance | none | **fixed** — per-call `gsl_root_fsolver_set_with_values` added (33b79b362); the `GSL_ERROR` handler contract is unchanged and documented (aaecb86e6) | gsl_root_fsolver_set "endpoints do not straddle y=0" |
| 63927 | 2023-03-14 | feature | Build | none | **fixed** — gsl-without-cblas.pc installed by both builds (442b7cbd2, 9ad4ed051) | Please add a gsl-without-cblas.pc |
| 64549 | 2023-08-12 | test | Performance | clean testcases.patch | **fixed** — one clean wrapper test replaces the four leaky near-duplicates (7b6979126) | diff patch containing new testcases for the interpolation module |
| 64613 | 2023-08-30 | - | Accuracy problem | inline | **fixed** — cdf beta_inc: 3 defects under fp-contract (1b1d94ee9) | fp-contract=fast stops convergence in beta_inc_AXPY/beta_cont_frac |
| 64777 | 2023-10-14 | bug | Runtime error | inline | rejected — not-a-bug: reporter's matrix is singular | gsl_linalg_complex_LU_decomp returns incorrect results |
| 64851 | 2023-11-03 | bug | - | clean inline gsl-config.in.patch | superseded — gsl-config exit status already done in the fork | gsl-config does not set correct status code |
| 65728 | 2024-05-12 | feature | - | none | **fixed** — complex `gsl_spblas_zgemv` and `gsl_spmatrix_*_set_add` (f218d8c83); complex real-vector scaling (3cdf34daa), sparse `_resize` (ae334370d), block/vector `_resize` (4a1b3d616) and complex `_unpack` (7e327561e) | Add sparse functionalities |
| 65760 | 2024-05-19 | bug | Accuracy problem | clean erfc.c.patch | **fixed** — erfc/log_erfc overflow rewritten (708791c25) | gsl_sf_log_erfc (and gsl_sf_erc) return NaN for infinite or very large finite arguments |
| 65868 | 2024-06-11 | - | - | clean 0001-bspline-Add-missing-definition-for-function.patch | **fixed** — gsl_bspline_eval_nonzero declared (ef3940132) | Missing definition for gsl_bspline_eval_nonzero |
| 65912 | 2024-06-23 | bug | - | clean 0001-Correct-GSL_SET_COMPLEX-if-native-complex-available.patch +1 | **fixed** — GSL_SET_COMPLEX signed zero / non-finite (631da98f8) | GSL_SET_COMPLEX is wrong if complex.h has no support for imaginary numbers |
| 65932 | 2024-06-30 | - | - | none | **fixed** - duplicate of #60457; the complex tridiagonal solvers are implemented (dd5e0a0f2) | Complex tridiagonal solvers - patch ignored? |
| 66026 | 2024-07-26 | - | Runtime error | partial bug_66026.patch +1 | **fixed** — LU_decomp_L3 identity-initialises ipiv (d73ba5300) | gsl_linalg_LU_decomp using uninitialized memory |
| 66128 | 2024-08-27 | doc | Documentation | clean specfunc.rst.patch | **fixed** — specfunc.rst result/val typo (066f5b747) | documentation of type gsl_sf_result: confused result and val |
| 66573 | 2024-12-18 | feature | - | clean vector_complex_conjugate.patch | **fixed** — `gsl_vector_complex_conjugate()` for the float, double and long double complex types, with a strided test (ec04a1312, d127fe3b2) | Feature: gsl_vector_complex_conjugate() |
| 66574 | 2024-12-18 | feature | - | clean 2_givens_cmplx.patch +3 | **fixed** — complex Householder right, Givens, bidiagonal decomposition and `gsl_linalg_complex_SV_decomp` (5df7c9876, b5325c1df, e0993e48c, ac4a51535) | Feature: complex singular value decomposition |
| 66575 | 2024-12-18 | feature | - | clean 5_svd_SV_solve_cmplx.patch +1 | **fixed** — `gsl_linalg_complex_SV_solve` and `_SV_decomp_mod` (77032d17c, 5a209fff9) | Feature: Extension to complex SVD |
| 66576 | 2024-12-18 | feature | - | clean fsolver_set_with_values.patch | **fixed** — `gsl_root_fsolver_set_with_values` for bisection/brent/falsepos, no ABI change (33b79b362) | Feature: gsl_root_fsolver_set_with_values() |
| 66695 | 2025-01-22 | feature | - | clean feagin_verner.patch | **fixed** — Verner 7(6)/8(7)/9(8) and Feagin 10(8)/12(10)/14(12) steppers (d67eae92b, 6fa6706f4) | Feature: additional high order ODE solvers (Feagin, Verner) |
| 66742 | 2025-01-31 | - | Documentation | clean gams.diff | **backlog** — large doc-only GAMS classification diff; not a correction, wanted but out of scope | GAMS classification |
| 66767 | 2025-02-08 | feature | Build | source binomialinv.c +1 | **fixed** — gsl_cdf_binomial_Pinv/Qinv by bisection, double result (e88a81c75) | Feature: inverse of binomial distribution |
| 66775 | 2025-02-10 | feature | Build | source poissoninv.c +1 | **fixed** — gsl_cdf_poisson_Pinv/Qinv by doubling + bisection, double result (e88a81c75) | Feature: inverse poisson distribution |
| 66800 | 2025-02-15 | feature | Accuracy problem | partial nbinomial-2.diff +1 | **fixed** — negative_binomial_pdf degenerate at p = 1 and p = 0 (957708351) | Feature: gsl_ran_negative_binomial_pdf with p = 1 |
| 66808 | 2025-02-17 | bug | Accuracy problem | source test_airy_zeroes_derivs.c | **backlog** — Airy accuracy not improved; sub-defects fixed (ac5f72d98); oscillatory accuracy beyond `eps*\|x\|^{3/2}` remains (residual defect below) | Bug: Airy Ai function values inaccurate |
| 66816 | 2025-02-19 | feature | Build | source nakagami_.c | **fixed** — gsl_ran_nakagami + pdf, sqrt of a gamma variate (999d53bbe) | Feature: Nakagami random distribution |
| 66826 | 2025-02-21 | feature | Build | none | **fixed** — test vectors added (728206d82) and the negative-integer-`b` defect they exposed fixed (75fcaf1b3) | Feature: test cases for function hyperg_1F1() |
| 66834 | 2025-02-23 | feature | Build | clean specfunc_gamma_test.diff | **fixed** — gamma_inc(0,0) = +Inf (062b15fc3) | Feature: gamma_inc(0, 0) handling and test cases |
| 66842 | 2025-02-24 | feature | Build | clean expint_infinity.diff | **fixed** — exponential integrals at the origin (9a00ba0ed) | Feature: exponential integrals at origin |
| 66844 | 2025-02-25 | feature | Build | none | **fixed** — dilog endpoint identities tested (25a16f928) | Feature: test cases for dilogarithm function |
| 66849 | 2025-02-26 | - | - | none | **fixed** — final step clamped and backed off; bsimp sub-step times computed directly (dda917296); RHS-records-the-time regression test | gsl_odeiv2_evolve_apply() may exceed final time |
| 66850 | 2025-02-26 | feature | Build | clean specfunc_hyperg_2F1.diff | **fixed** — 2F1 at c = a-1 and c = b-1 with the corrected relation (adc9389a4); the posted patch is mathematically wrong | Feature: special cases hypergeometric2F1(a+1, b, a, x) and hypergeometric2F1(a, b+1, b, x) |
| 66862 | 2025-03-02 | bug | - | partial specfunc_test_sf.diff +1 | **backlog** — report not a bug (duplicate complex sin/cos; submitter withdrew); sharing the complex sin/cos implementation between `complex/math.c` and `specfunc/trig.c` is wanted refactoring | Bug: duplicate functions in GSL produce differing results |
| 66874 | 2025-03-05 | feature | - | source doc_bst.rst | rejected — superseded: the bst module is documented in the fork (20eda0316) | Feature: documentation for binary search trees (bst module) |
| 66877 | 2025-03-06 | feature | Build | none | **fixed** — exact zero error at x = 0 (fe7285909) and the special-case vectors (2d6cf433c) | Feature: test cases for hyperg_0F1() |
| 66880 | 2025-03-07 | feature | Build | clean test_hyperg_U.diff | **fixed** — U special-case vectors, with corrections (153a149d4) | Feature: test cases for function hyperg_0F1() |
| 66886 | 2025-03-09 | feature | Build | clean specfunc_trig_.diff | **backlog** — refactor of `specfunc/trig.c`, not a bug fix; wanted but out of scope | Feature: refactoring of specfunc/trig.c |
| 66894 | 2025-03-11 | - | Documentation | none | **fixed** — chapter introduction notes that randist does not validate distribution parameters (dbe316d49) | Domain value checking for random number distributions |
| 66922 | 2025-03-17 | feature | Build | clean specfunc_trig.diff | **fixed** — complex_logsin at a real zero of sin (10d64fb15) | Feature: special case for special functin logsin() |
| 66949 | 2025-03-25 | feature | Build | source erlang.c | **fixed** — gsl_cdf_erlang_P/Q delegating to the gamma cdf (3d607fce8) | Feature: Erlang cumulative distribution |
| 66993 | 2025-04-05 | bug | Accuracy problem | none | **backlog** — zero-base overflow documented (afd7f7d6c); the other claims do not reproduce (effectively closed) | Bug: pow_int issues |
| 67058 | 2025-04-28 | - | Accuracy problem | none | **fixed** — empty data sets raise `GSL_EBADLEN` (NaN with the handler off) (68fc3c752) | gsl_stats_mean and gsl_stats_sd result to 0.0 when array with zero length |
| 67301 | 2025-07-10 | bug | - | none | rejected — not-a-bug: r2 is a within-bin fraction, a second find() would pick the wrong bin; documented 2b15efbff | Bug: Test for existence of uniform random variate in histogram |
| 67359 | 2025-07-27 | feature | - | clean hyperg_2F0.c.patch +1 | **fixed** — 2F0 with a zero parameter (1d99544d9) | Feature: Add special case for gsl_sf_hyperg_2F0 |
| 67445 | 2025-08-20 | - | Build | clean 0001-linalg-increase-cholesky_invert-Hilbert-test-toleran.patch | **fixed** — cholesky_invert Hilbert tolerance (f84a57a0f) | linalg test failures under gcc 14.2.1 |
| 67446 | 2025-08-20 | - | Build | none | **backlog** — gcc 14.2.1; passes on gcc 15.2 x64; re-check parked pending gcc 14.2.1 | multilarge_nlinear test failures under gcc 14.2.1 |
| 67447 | 2025-08-20 | - | Build | none | **backlog** — gcc 14.2.1; passes on gcc 15.2 x64; re-check parked pending gcc 14.2.1 | spmatrix test failures under gcc 14.2.1 |
| 67494 | 2025-09-09 | - | Accuracy problem | none | rejected — cosmetic; `gsl_pow_int` and `gsl_sf_pow_int` return identical values | Consistency in pow_int usage |
| 67621 | 2025-10-23 | - | Documentation | clean rst_deprecated.diff | **fixed** — deprecation markers from headers (5df6c0079) | Mark functions as deprecated if they are also deprecated in the code |
| 67689 | 2025-11-10 | doc | Documentation | clean specfunc-psi.rst.patch | **fixed** — gsl_sf_complex_psi_e documented (345883172) | Incomplete documentation of digamma functions in GSL specfunc |
| 67705 | 2025-11-15 | bug | Build | none | rejected — not reproducible; linalg passes with a fixed seed | ttest failure in linalg/QR_solve_r random |
| 67728 | 2025-11-23 | bug | Accuracy problem | source psi_dropin.c | **fixed** — polygamma at negative arguments (e9f69933a) | gsl_sf_psi_n_e yields domain error |
| 67774 | 2025-12-05 | bug | Accuracy problem | none | **fixed** — the value was already right for x < 0; the error bar was computed with the signed x and came back negative (4220f05d9) | Feature: arctan integral is also defined for negative inputs |
| 68068 | 2026-02-19 | bug | Accuracy problem | none | **fixed** - quad_golden stored f_m in f_upper instead of f_lower (4260778fe) | Bug: incorrect straddling of area of convergence in quad_golden |
| 68073 | 2026-02-20 | bug | Documentation | none | **fixed** — gsl_stats_select comment corrected from "k-th largest" to "k-th smallest" (521bed61d) | Bug: incorrect inline code comment on BASE FUNCTION(gsl_stats,select) |
| 68098 | 2026-02-27 | feature | - | none | **fixed** — \|rho\| >= 1 rejected with GSL_EDOM in the bivariate Gaussian pdf and generator (6da0c08f3) | Feature: division by zero when data is perfectly correlated |
| 68283 | 2026-04-26 | - | - | none | **backlog** — the proposed guards break the `gsl_stats` equivalence in `rstat_test`; a consistent skew/kurtosis correction remains | Correction to gsl_rstat_skew and gsl_rstat_kurtosis |
| 68312 | 2026-05-07 | bug | Accuracy problem | none | **fixed** — same recurrence as #43256 (59fc479e2) | Wigner symbols inaccurate for large j |
| 68367 | 2026-05-19 | feature | - | source invelljac.c | **fixed** — `gsl_sf_elljac_arcsn_e/arccn_e/arcdn_e`, using Carlson RF with the parameter `m`; the posted prototype passed `m` where the modulus `k = sqrt(m)` was expected (183ca894f, 653c77bca) | Feature: inverse Jacobi elliptic integrals |
| 68379 | 2026-05-21 | - | - | none | **fixed** — size_t for the accessor index (ae60e2883) | Histogram: expand scope of internal variables |
| 68398 | 2026-05-26 | feature | - | none | **fixed** - eta_int leading term via gsl_ldexp, value correctly rounded, error halved (f774247e7) | Feature: use GSL native gsl_ldexp in eta fuction for integer argument |
| 68415 | 2026-06-02 | bug | - | none | rejected — not-a-bug: variance/covariance/pvariance all return NaN at n=1; n >= 2 documented (65bc71664) | Bug: inconsistency in variance error handling in statictics module |
| 68479 | 2026-06-25 | bug | Accuracy problem | none | **fixed** — BINV seed as `exp(n log1p(-p))` (7376dc316) | Bug: gsl_ran_binomial is not accurate |
| 68495 | 2026-07-03 | bug | Build | inline | **fixed** — pow_int INT_MIN overflow (41b1e2c00) | Undefined behavior in gsl_pow_int: signed int overflow |
| 68518 | 2026-07-13 | - | Build | none | **fixed** — DISCARD_POINTER replaced by void casts and the macro dropped from both configs (46fa69ae1, 2afdbb721) | Remove stale construct in configure.ac |
| 68549 | 2026-07-23 | - | - | none | **backlog** — uniform argument lists across the distribution functions; breaking API change, wanted but out of scope | Breaking change: match argument list among distribution functions |
| 68592 | 2026-08-03 | doc | Documentation | none | **fixed** — new BST chapter added to both manuals and to the index/menu, marked as an unreviewed AI draft (20eda0316) | Documentation: missing Binary search tree documentation from index |
| 68611 | 2026-08-13 | feature | - | none | **fixed** - all-empty histogram raises GSL_EDOM in pdf_init (b60a47106) | Feature: guard against degenerate histogram input |
| 68625 | 2026-08-21 | bug | Accuracy problem | none | **fixed** — hermite_func_der at n=0,1 (5d36e1999) | Bug: incorrect derivative of the Hermite function of order 0 or 1 |
| 68663 | 2026-08-31 | feature | - | none | **backlog** — drop the `GAMMA_INC_A_0` alias (one line); no correctness case, low value | Feature: remove GAMM_INC macro in favor of direct call |
| 68704 | 2026-09-18 | feature | - | none | **backlog** — `gsl_histogram_variance` (new API) including the breaking rename question | Feature: introduce variance of a histogram gsl_histogram_variance |


## Future work backlog

Every open item has a verdict, and the eligibility rule has since been
widened to **any report on which an action can be taken**.  The only
reports left out are those about platforms this fork cannot build and test
(32-bit targets and architectures absent from CI), plus reports that name
no defect in the repository at all.  The backlog therefore carries, besides
the feature/API tail, the residual work on reports whose reported defect is
only partly fixed or declined, and several open defects that until now were
recorded only in `SAVANNAH_REVIEW.md`.  These rows carry the **backlog**
verdict in the table above, kept distinct from **rejected**.

Nothing here is a fork commitment; items are removed one by one as they are
taken.  This file is scratch and is not part of the fork's record of
changes, so nothing here needs a `FORKNEWS` entry until an item is actually
implemented.

### History of the widening

The first stage admitted only bug fixes, documentation corrections and
test-quality improvements.  Each later widening is recorded below so the
reason an item sits in the backlog is not lost.

Four performance-only items that were listed here have since been taken as
result-preserving fixes once the filter was widened: `#54925` (cblas gemm
loop reorder), `#51104` and `#21833` (`gsl_permutation_next`) and `#40092`
(false-position redundant evaluation).  See `FORKNEWS`; they are no longer
listed in the table.

Five test-quality items have also been taken: `#66826`, `#66844`, `#66877`,
`#66880` and `#64549`.  See `FORKNEWS`; they are no longer listed in the
table.  One of them, `#66826`, exposed a separate library defect in
`gsl_sf_hyperg_1F1_int_e` for negative integer `b`, recorded as an open item
in `SAVANNAH_REVIEW.md`; the fix is not part of the test change.

Six items previously rejected here as "new special cases" were taken once
the eligibility filter was widened to admit correctness fixes to existing
API that add no new symbol: `#66800`, `#66834`, `#66842`, `#66850`,
`#66922` and `#67359`.  Each is a defect for inputs the library already
claimed to handle, so all six are now fixed; they are no longer listed in
the table.  See `FORKNEWS`.

Three further deferred items were taken once the filter was widened again:
`#30324` (2F1 continued to `x < -1`), `#30540` (convergence checks in the
v1 `rk2imp`/`rk4imp`) and `#68098` (the bivariate Gaussian `rho` range).
They are no longer listed in the table.  See `FORKNEWS`.

The filter was widened once more to admit new public API.  Five items
listed below were taken: `#59900` (truncated Gaussian), `#66816`
(Nakagami), `#66949` (Erlang CDF), `#66767` (inverse binomial) and
`#66775` (inverse Poisson); `#24252` (gamma tail) was re-examined and
**rejected** because the posted sampler is numerically wrong (see
`FORKNEWS`).  `#24871`, the large new complex algorithm, was taken
later as a port of Amos Algorithm 683.  All six are no longer listed in
the table.

Two further new-API items have since been taken: `#66573`
(`gsl_vector_complex_conjugate()`) and `#68367` (inverse Jacobi
elliptic functions, where the posted prototype misused the modulus).
`#60457` and its duplicate `#65932` (complex tridiagonal solvers) were
taken as well, followed by the complex SVD suite `#66574` and `#66575`,
the Feagin/Verner ODE steppers `#66695`, the configurable
finite-difference Jacobian step `#45782`, the complex zeta/eta/Hurwitz
functions `#57173` and the new fixed quadrature rules `#32257`.  They
are no longer listed in the table; see `FORKNEWS`.

`#41527` was taken as a bug fix once the filter was widened to admit
contract changes to existing API: the multimin simplex minimizers now
return `GSL_EBADFUNC` for a non-finite objective instead of masking it as
`GSL_EFAILED` (`f007aa42c`); it is no longer listed in the table.

`#31109` was taken as a result-preserving-enough performance fix once the
`ode-initval2` framework's driver/control link was used: `bsimp` derives
its extrapolation order from the driver control's error level, so a loose
tolerance no longer pays for the highest order (`1fbb40ca3`).  It is no
longer listed in the table.

`#65728`, already partly fixed (`f218d8c83`), is now complete: complex
real-vector scaling, `gsl_spmatrix_resize`, `gsl_block_resize`/
`gsl_vector_resize` and `gsl_vector_complex_unpack` were added
(`3cdf34daa`, `ae334370d`, `4a1b3d616`, `7e327561e`).  It is no longer
listed in the table.

A defect in the same 1F1 family was found while fixing `#28267`: for a
non-positive integer `a` with `b > 0` and large `|x|`, the terminating
value went through two nearly reciprocal Kummer factors, which can
under/overflow while the product is representable, so
`gsl_sf_hyperg_1F1_int_e` returned `0` or `NaN` with `GSL_SUCCESS`
(`1F1(-100, 2, 1000)` against a true `1.2808e135`), and the
`x -> +Inf` asymptotic was reached for a negative integer `a`, where it
is invalid.  It has no Savannah report; fixed in `43dedeaed` (see
`FORKNEWS`).

`kind`: `feat` = new feature/API, `perf` = performance only, `test` =
test-quality request, `doc` = documentation matter, `build` = packaging or
distribution, `refactor` = internal restructuring, `accuracy` = accuracy
improvement, `assess` = needs a fresh verdict.  `patch`: `clean`/`partial`/
`source` if an attachment was collected by the sweep, `-` if none;
`inline` if the text is in the bug thread only.

### Residual defects

Real defects, or work left open by a partly applied fix.  Most are already
analysed in `SAVANNAH_REVIEW.md`; the reader is referred there for the
reproduction.  `#66826` had no triage row of its own: the report was closed
test-only, but the library defect it exposed was later fixed (75fcaf1b3).

| # | kind | what remains |
|---|---|---|
| 21835 | bug | 2F1 at `c = a+b` near `x = 1` exhausts the series (`GSL_EMAXITER`) despite ~9 good digits; needs the A&S 15.3.10 form. |
| 50712 | bug | `lm+accel` with a finite-difference `fvv` on `box3d` fails to converge; reproduces on MSVC x64 and is disabled under `#if 0`. |
| 43259 | assess | re-verdict the `pochrel` and `ellint_P` inputs dismissed as "degenerate"; they may be a genuine limitation rather than evaluation on a pole. |
| 66808 | accuracy | Airy oscillatory accuracy beyond `eps*\|x\|^{3/2}` (a double-double prototype recovered only 5-10x, capped by the double coefficient tables). |
| 36152 | test | `gsl_sf_sinc_e` large-`x` vector (the behaviour is fixed as a side effect; no vector yet). |
| 66993 | doc | pow_int residual, effectively closed: the only valid point, the zero-base overflow, is documented; the other claims do not reproduce. |

### Wanted work, no defect

Feature, API, performance, refactor and documentation items the earlier
rule set aside and the widened rule admits.

| # | kind | patch | what it would add |
|---|---|---|---|
| 25320 | feat | - | Fresnel integrals, with the extension's negative-`x` sign defect corrected. |
| 21837 | feat | - | permutation path for a zero-diagonal symmetric tridiagonal solve (new algorithm). |
| 34361 | feat | - | `gsl_bspline_knots_greville` constrained least squares (Lawson-Hanson NNLS/LDP). |
| 68704 | feat | - | `gsl_histogram_variance`, including the breaking rename question. |
| 24252 | feat | source | a *correct* truncated-gamma / gamma-tail generator (the posted sampler is wrong). |
| 50382 | build | - | NuGet packaging for Windows. |
| 59914 | doc | - | a documented native Windows build write-up. |
| 66862 | refactor | partial | share the complex sin/cos implementation between `complex/math.c` and `specfunc/trig.c`. |
| 68663 | refactor | - | drop the `GAMMA_INC_A_0` alias (one line; low value). |
| 29834 | refactor | source | unify the CBLAS and GSL argument-checking macros. |
| 68283 | bug | - | an `gsl_rstat` skew/kurtosis correction consistent with the `gsl_stats` equivalence the module asserts. |
| 68549 | feat | - | uniform argument lists across the distribution functions (breaking). |
| 47402 | perf | - | Mathieu coefficient caching in the workspace (items 1, 2 and 4 already done). |
| 66742 | doc | clean | GAMS classification across 177 files. |
| 66886 | doc | clean | refactor of `specfunc/trig.c`. |

### Platform re-checks

Reports that named a platform or toolchain the main review did not exercise.
The first two are reachable from the current CI/setup; the compiler-specific
ones need that toolchain installed first.

| # | what to do |
|---|---|
| 56843 | reproduce the non-x86 eigen accuracy difference on the arm64 macOS CI target, which is in CI. |
| 44865 | make the `e5_bigt` comparison absolute-aware so it survives `-ffp-contract=fast`; the contraction path itself is testable on this machine. |
| 67446, 67447 | rebuild with gcc 14.2.1 and attempt reproduction. |
| 49697 | re-verify with `-mavx`; likely already covered by the raised cholesky tolerance (`#67445`). |
| 39152, 54919 | attempt reproduction with Intel oneAPI `icc`/`icx`. |

### Excluded

Skipped because the platform cannot be built or tested here: 32-bit targets
(`46593`, `52322`, `59759`, `53903`, plus the 32-bit glibc data points
`45265` and `52927`) and architectures absent from CI (`47028` ppc64le,
`48915` AIX, `61342` ppc64/sparc).  `42058` is Savannah's own feed and
`45797` is an external LAPACK build; neither is a defect in this
repository.  `66874` is superseded - the fork already documents the BST
module in `20eda0316`.

The rest of the main table is handled or needs no action; see the per-bug
verdicts and `FORKNEWS`.
