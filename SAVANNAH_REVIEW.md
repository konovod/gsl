# Savannah bug review — working notes

Scratch notes for the review of the patches published in the Savannah GSL
bug tracker (https://savannah.gnu.org/bugs/?group=gsl).  This file records
the review verdicts so the work is not lost; it is **not** part of the
fork's record of changes and is deliberately kept out of `FORKNEWS`,
which only carries changes that were actually made.

Applied changes and their reasoning live in `FORKNEWS`.


## How the inventory was produced

`scripts/savannah_bugs.py` walks the open bug list, then each bug page,
downloads the attached patch-like files, and runs `git apply --check`
against the working tree.  No login is needed; attachments come from
`file.savannah.gnu.org`.

    python scripts/savannah_bugs.py            # writes to %TEMP%/sav-gsl

A full sweep is ~220 requests and takes about eight minutes.  Results are
cached in `<out>/cache.json`, so re-running the report is instant.

    %TEMP%\sav-gsl\inventory.tsv      # 219 rows, 11 columns
    %TEMP%\sav-gsl\patches\<bugid>\   # downloaded patches
    %TEMP%\sav-gsl\cache.json         # parsed bug pages

Counts from the 2026-10-01 sweep:

| | count |
|---|---|
| open bugs | 219 |
| carrying a patch (attached or inline) | 88 |
| ... of those, already handled by this fork | 8 |
| `git apply` clean | 42 |
| partially applies (several attachments, mixed) | 8 |
| inline only, nothing to download | 9 |
| whole replacement files, not diffs | 28 |
| does not apply | 34 |

Note that "most recent 10 bugs" is a bad window: only one of the newest
fifteen carries a patch, against roughly two thirds of bugs 16 to 79.
The recent submissions are mostly feature requests with the discussion
continued on the mailing list.


## The main lesson

`git apply --check` only proves a diff still applies.  It says nothing
about whether the change is *correct*, and several Savannah patches
apply perfectly while being wrong or actively harmful.  Every candidate
has to be tested against the built library.  See `#52321` below.


## Applied

| bug | verdict |
|---|---|
| `#65868` | applied as posted (commit ef3940132) |
| `#43326` | applied, one hunk dropped as dead code (commit 4f9f4f4fc) |
| `#67445` | applied as posted (commit f84a57a0f) |
| `#66128` | applied as posted (commit 066f5b747) |
| `#57978` | applied with `GSL_EDOM` instead of the patch's `GSL_ELOSS`, and the missing documentation added (commits 7c27b358b, f87e96143, b4a012ca9) |
| `#53919` + `#65760` | applied together, rewritten: neither patch fixes the actual overflow (commits 708791c25, 6bc4d8cdd, 9daf0cf36) |
| `#66862` | rejected as not-a-bug; the test-coverage gap it exposed was filled separately (commits 48a49003f, 758c15422) |
| `#64613` | reproduced on this machine (it has FMA) and fixed: three separate defects, not the one reported (commits 1b1d94ee9, 98e966458) |
| `#39292` | **fixed 2026-10-02** (commit 90d9037a5): the rejection did not reproduce and an independent check favours the patch; `C = 0.5 sqrt(1 + 4Q)` is applied and a turning-point guard vector added.  Earlier verdict, for the record: rejected — the posted `0.5` factor makes the branch *worse* (commit 192ec9cbf) |
| `#64777` | rejected - the reporter's matrix is singular and the permutation was ignored |
| `#68625` | fixed - `gsl_sf_hermite_func_der_e` dropped the `-x psi_n` term for `m=1, n=0,1`, and reported `err = 0` at `n = 0` (commits 5d36e1999, f4d35a67e).  Follow-up: the exact-zero `psi_1'(1)` vector failed on arm64 at `TEST_SF_INCONS`; the `m = 1` error estimate was fixed to include the cancelling terms and the vector now passes at `TEST_TOL1` (commits 491f38efa, a71516e66) |
| `#67689` | applied - `gsl_sf_complex_psi_e` was undocumented; patch taken with three corrections (commit 345883172) |
| `#58066` | **mis-assessed first time round**; re-examined, the reporter was right, and the fix is applied to `psi`, `psi_1` and `complex_psi` (commit d53d00109) |
| `#47345` | applied - `gsl_complex_arccosh` returned `-0` for a real argument; now delegates to `gsl_complex_arccosh_real` (commit d6ec47d87) |
| `#60371` | applied - the 2D interpolation wrapper left the output unwritten on a domain error; now stores `NaN` (commit 2b4e2f1e5) |
| `#54077` | applied - `gsl-randist` read its seed with `atol`, saturating large seeds; now `strtoul` (commit aebe57a5c) |
| `#66026` | applied with a different fix than posted - `LU_decomp_L3` leaves `ipiv` unwritten on a singular pivot; the array is identity-initialised rather than skipping the permutation (commit d73ba5300) |
| `#52359` | applied - the Airy modulus/phase error estimates divided by series corrections that can vanish; rewritten from the full expressions (commit ac5f72d98) |
| `#52570` + `#51000` | applied - huge negative arguments gave `±inf`/`NaN` or non-physical magnitudes with status `GSL_SUCCESS`; now `GSL_ELOSS` with a bounded result (commit ac5f72d98) |
| `#66808` | partly applied - the oscillatory accuracy is phase-limited at `~eps*|x|^{3/2}`; a double-double phase was prototyped and rejected (the double-precision coefficient tables cap it), and the vectors are added at realistic tolerances (commit ac5f72d98) |
| `#58067` | not-a-bug - `Ai(113)` genuinely underflows; `GSL_EUNDRFLW` is correct and is now documented (commit 86cb80782) |
| (Olver) | the Airy fix exposed an optimistic error estimate in `bessel_olver.c`; the missing argument-error term is now propagated (commit 790fe1078) |
| `#67621` | applied with corrections - the posted patch touched only two files and marked the wrong Hermite functions; the markers are derived from the actual `GSL_DISABLE_DEPRECATED` blocks.  `gsl_bspline_knots_greville` is deliberately **not** marked: it is only commented "future to be deprecated", not guarded (commit 5df6c0079) |
| `#43902` | applied - `doc/examples/vectorw.c` now writes 10 elements, matching `vectorr.c` (and its `doc_texinfo` mirror) (commit 0e82d8223) |
| `#44952` | applied - the unreachable `<varargs.h>` branch in `test/results.c` is folded into the `<stdarg.h>` branch (commit e2354de55) |
| `#42472` | applied with a larger fix than posted - `HH_svx`/`HH_solve` had the row/column test backwards and `HH_solve` overflowed `x` for `b` longer than `x`; both now share a tall Householder QR (commit 1f84bed77) |
| `#65912` | **fixed 2026-10-02** (commit 631da98f8): the native branch built `x + I*y`, which loses the sign of a zero component and turns a non-finite part into `NaN`; the macro now assigns the components directly through `GSL_REAL`/`GSL_IMAG`.  Earlier verdict, for the record: rejected - both posted variants (`CMPLX`, `_Generic`+`CMPLXF/L`) are unusable |
| `#59834` | **fixed 2026-10-02** (commit 4e4a88242): the posted patch does actually align `maxque`, contrary to the earlier note; the fork instead adds an explicit `ringbuf_align()` helper and fixes the same defect in `qnacc` found while auditing |
| `#47646` | code fix already upstream (`05c5b5179`); this fork only adds the missing regression test (commit d64cc4d93) |
| `#36152` | **fixed 2026-10-02** (commits a08ef2f7f, 8a46ec7cf): the report is about the spherical Bessel family, and the *j* half was already fixed upstream (`cd2dd0519`, `bd5b94b47`).  The surviving defect was (a) the Y functions still calling `gsl_sf_sin_e`/`cos_e`, and (b) the underlying reduction in `gsl_sf_sin_e`/`cos_e` itself; both are fixed.  Also closes `#45726` and the trigonometric half of `#45746` |
| `#68495` | **fixed 2026-10-02** (commit 41b1e2c00): the reported `-n` at `gsl_pow_int` is UB for `n = INT_MIN`; fixed as posted, and the identical negation in `gsl_sf_pow_int_e` - which `9493ac014` missed and which loops forever - is fixed with it |
| `#32306` | **fixed 2026-10-02** (commits e4c4ac326, 882c8361d): the integer-`c-a-b` branch of `hyperg_2F1_reflect` forms every gamma factor with `gsl_sf_lngamma_e()` and applies a single global sign, so `2F1(-1/2,3/2;1;x)` came back negated for `x >= 1/2`; integer-`d` cases with `x < 0.995` now use the Gauss series instead.  The `err = 1` for a one-signed series (`a < 0`) is fixed as well.  Also covers #54998 and the Monajemi case of #39056.  Residual: reflection at `x >= 0.995` keeps its sign/accuracy defects |
| `#43809` | **fixed 2026-10-02** (commit 52505315d): the `a < 0, b > 0` branch reduced to Kummer and evaluated the transformed call with an unstable backward recurrence on `b`; for the reported parameters that returned `3.39e80` instead of `7.51e60`.  The direct series is now evaluated too and preferred when it reports `err / abs(val) < 1e-10`.  Regression test added. |
| `#28267` | partly addressed by the `#43809` fix — the transition region `x ~ abs(a)^2` still loses most digits (e.g. `(-37.8, 2.01, 103.58)` at ~2%); neither the recurrence nor the series is accurate enough there, so it stays open |
| `#39372` | **fixed 2026-10-02** (commit c1df353ae): `gsl_hypot3` divided by `max(|x|,|y|,|z|)`, so an infinite argument produced `inf/inf = NaN`; the function now returns `+Inf` for any infinite argument. The reporter's suggested `gsl_hypot(gsl_hypot(x,y),z)` was not used: `gsl_hypot` raised a range error for infinities at the time (the companion `#57979` fix) |
| `#57979` | **fixed 2026-10-02** (commit 1a470222a): `gsl_sf_hypot(NaN, y)` returned `sqrt(2)*y` and `gsl_sf_hypot(inf, y)` reported overflow, because `GSL_MIN_DBL`/`GSL_MAX_DBL` never select a NaN; both cases now follow the C99 spec (`+Inf` / NaN, `GSL_SUCCESS`) |
| `#37894` | already fixed upstream, inherited unchanged — commits `ae19e3e8b` (`LT_INIT([win32-dll])`, `LT_LIB_M`) and `1d002ee93` ("add cygwin patch from J.P. Flori") turn the old MinGW-only conditionals into a `*-*-cygwin* | *-*-mingw*` test that adds `-no-undefined` and `GSL_LIBADD=cblas/libgslcblas.la` to the shared libraries.  Both commits are ancestors of `savannah/master`, so the posted `gsl-autotools.diff` is fully present; no fork change is made |
| `#58032` | **fixed 2026-10-02** (commit 1f2607a0b): `gsl_sf_hyperg_1F1_e`/`_int_e` tested `x == 0` before any parameter check and returned 1, so the pole at `b = 0, -1, -2, ...` was reported as `1F1(a,b,0) = 1`.  The branch now classifies the parameters as the `x != 0` branches do and returns `GSL_EDOM` there, while terminating cases (`a` a nonpositive integer with `a >= b`) stay at 1.  Manual updated in both trees.  No patch was posted to the tracker, so this is a from-scratch reproduction |
| `#39057` | **fixed 2026-10-02** (commit 9befdae95): the reported call is real but its expected value, `0.99477710813146`, is `gsl_cdf_chisq_P(0.5, 0.01)`, not the inverse.  The gamma inverse's branch heuristic is only valid for `a >= 1`; with `a = 0.005` the fixed-step iteration cannot reach `3.51e-61`, and `gsl_cdf_gamma_Qinv` returned wrong values silently.  The solver is reworked around a bracketed Pegasus iteration on `log x`, and both entry points delegate to the smaller tail.  No patch was posted |
| `#53451` | **fixed 2026-10-02** (commit 2673ce0d8): the title says `Pcomp` but the report's numbers are the *incomplete* `gsl_sf_ellint_P`.  Both it and `Pcomp` evaluate `Pi` through `RJ(..., 1 + n sin^2(phi))`, which rejected `p < 0` with `GSL_EDOM`; that is exactly the case where the integrand has a pole and the value is the Cauchy principal value.  `gsl_sf_ellint_RJ_e` now applies the DLMF 19.20.14 / Boost.Math transformation for `p < 0` (`p = 0` still `GSL_EDOM`), and `gsl_sf_ellint_P_e`'s error estimate no longer propagates a signed `n/3`.  Five vectors added; manual updated in both trees |
| `#43256` + `#68312` | **fixed 2026-10-03** (commit 59fc479e2): the direct Racah 6j sum overflowed at about `171!` (the report's `(14,16,16;78,62,76)` is the doubled form of `(28,32,32;156,124,152)`) and, below that, cancelled catastrophically - at `j = 100` the largest term is `~3e14` times the result.  Replaced by the Schulten-Gordon recurrence (public-domain SLATEC `DRC6J`), with backward/forward matching and orthogonal normalisation; 9j and Racah W inherit it.  References from sympy and an exact Racah evaluation; negative control fails with `val = inf`, `GSL_EOVRFLW` |
| `#45924` | **fixed 2026-10-03** (commit 51a63cbd5): the beta inverse bisected to an absolute `Ptol = 0.01` and then ran an unsafeguarded Newton iteration, so large `a,b` at any tail away from one half (and `a << 1`, and `Qinv`'s `1 - Pinv` complement) returned `NaN`.  Reworked on `t = logit(x)` with a maintained bracket, safeguarded Newton and direct upper-tail evaluation.  Verified against scipy `betaincinv`/`betainccinv` over 173000 points, worst relative error `7.4e-11`.  Ten report vectors plus fdist delegation vectors added |
| `#67058` | **fixed 2026-10-03** (commit 68fc3c752): `gsl_stats_mean` returned 0.0 for an empty data set (the recurrence leaves the accumulator at zero) and the variance/sd family returned 0.0 through it; R and numpy return `NaN`.  `gsl_stats_mean` and `compute_variance` now raise `GSL_EBADLEN`, returning `NaN` with the handler off; `tss` is left at 0.  Both manuals updated.  99 new failures with the guards reverted |
| `#42502` | **rejected** (not a bug): `gsl_cdf_ugaussian_Pinv(0.5)` returns exactly 0.0 on the built DLL and is covered by `cdf/test.c:498`.  The reporter's program omits `<gsl/gsl_cdf.h>`, so the function is implicitly declared `int` and `printf`'s second `%f` reads a stale vararg slot (the `1.000000`) |
| `#66844` | applied - two exact dilog endpoint identities added at `TEST_TOL0` (commit 25a16f928) |
| `#66826` | test vectors applied (commit 728206d82); the negative-integer-`b` defect the report exposes is **open**, see the section below |
| `#66877` | applied - `gsl_sf_hyperg_0F1_e` now reports an exact zero error at `x = 0` (commit fe7285909), and the `x = 0` / half-integer-`b` vectors are added (commit 2d6cf433c) |
| `#66880` | applied with corrections - `M_CUBEROOT2` written out, two vectors the reporter disabled are enabled, two inaccurate ones left out (commit 153a149d4) |
| `#64549` | applied as a single clean wrapper test instead of the four leaky near-duplicates (commit 7b6979126) |

## Resolved

Note: `#57978` and `#53919`/`#65760` resolve in opposite directions on
purpose.  `sin_pi` has no limit at infinity, so an infinite argument is
out of domain (`GSL_EDOM`); `erfc` and friends all have finite limits, so
infinity is in the domain (`GSL_SUCCESS`).  Getting that backwards in
either direction would have been wrong.

### `#53919` + `#65760` - erfc / log_erfc / erf_Z at huge arguments - applied

Both bugs describe "NaN for infinite or very large finite arguments", and
the reported symptoms are real, but **my first reading of the cause was
wrong** and the value — not only the error estimate — was corrupted.
Recording that because the patches both look plausible:

* I assumed the values were already correct and only `err` was NaN,
  reasoning that `erfc8_sum()` "reproduces the correct asymptote
  indefinitely". It does not. It is evaluated by Horner to degree 5 and
  degree 6 (`erfc.c:69-76`), so **both** `num` and `den` overflow to
  infinity above `x ~ 1e62` and the ratio is `NaN`. That is what actually
  produced the `NaN` value.
* I then "fixed" only the error term, rebuilt, and re-ran: `erfc`,
  `log_erfc` and `erf` were *still* `NaN` for large `x`. The guard I had
  written (`e_val == 0.0`) never fired.
* Second, subtler gap: my `log_erfc` fallback was keyed on
  `gsl_isnan(e_val)`, but in the window `2.4e51 < x < 1e62` `den` has
  overflowed while `num` has not, so the ratio underflows to `0` and
  `log_erfc8` returns `-inf`, *not* `NaN` — while the true answer
  (`-5.8e102` at `x = 2.4e51`) is perfectly representable. Needed
  `gsl_finite()`, not `gsl_isnan()`. Caught only by sweeping a table of
  values across the transition instead of testing the endpoints.

Where the real limits are, all measured:

| x | mechanism |
|---|---|
| `27.213` | `erfc(x)` reaches the least denormal, 4.9e-324 |
| `2.4e51` | `erfc8_sum()` denominator overflows; ratio -> 0, `log` -> `-inf` |
| `1e62` | numerator overflows too; ratio -> `NaN` |
| `1.34e154` | `x*x` overflows, so `-inf` is the honest answer for `log_erfc` |

Final behaviour, verified by running the built library:

    erfc_e   ( inf)      -> val=0          err=0
    erfc_e   (-inf)      -> val=2          err=8.8818e-16
    log_erfc_e( inf)     -> val=-inf       err=0
    log_erfc_e(-inf)     -> val=0.6931472  err=7.5191e-16
    erf_e    ( inf)      -> val=1          err=4.4409e-16
    erf_Z_e  ( inf)      -> val=0          err=0   (GSL_EUNDRFLW)
    log_erfc_e(2.4e51)   -> val=-5.76e+102 err=2.558e+87
    log_erfc_e(1e52)     -> val=-1e+104    err=4.441e+88
    log_erfc_e(1.5e154)  -> val=-inf       err=0

Seamlessness of the switch to the asymptotic form was checked by
comparing against `-x*x - log(x) - LogRootPi_` over a sweep: relative
difference is exactly 0 from `2e51` through `1.3e154`, no discontinuity.

Why neither patch was usable as posted:

* `#53919` (`erf.patch` and `v2-erf.diff`, near-identical) sniffs for
  `isnan` in the *result* and guards with `!isnan(x)`. It never addresses
  the polynomial overflow, so it only converts `NaN` into a wrong finite
  answer for the `x` window where the ratio underflows to `0`. Its
  `else`-branch comment "negative overflow" is also mislabelled — `erfc`
  underflows there, it does not overflow.
* `#65760` has a good `erfc` half (its `#if 0` block documenting the
  `27.213` / `26.213` underflow points is the right observation) but a
  bad `log_erfc` half: it discards `log_erfc8()`, which is still
  accurate well past its chosen threshold of `x = 100`, replacing it with
  `-x*x - log(x) - M_LNPI/2` — the same asymptotic term, since
  `M_LNPI/2 == LogRootPi_` — with a hard-coded `err = 0.0` over an
  unbounded range. It also misses the `x > 1.37e154` overflow its own
  comment describes.

Two decisions taken, both confirmed by the user:

* **Include `gsl_sf_erf_Z_e`**, which neither patch mentions. It had the
  same `inf * 0` shape in `fabs(x * val)` at `erfc.c:382`. Worth noting
  it was broken at `+-inf` only, not for large finite `x`, since
  `exp(-x*x/2)` underflows before `x * val` can overflow.
* **`GSL_SUCCESS`**, not `GSL_EOVRFLW`, for `log_erfc` beyond
  `sqrt(DBL_MAX)`. `doc/specfunc-erf.rst` said `Exceptional Return
  Values: none`, so inventing an exceptional return would have
  contradicted the documented contract; the docs were corrected instead.

Rejected both patches' `!isnan(x)` clauses, consistent with the `#57978`
decision: no `specfunc` function treats a `NaN` argument as a domain
error, it propagates.

### `#64613` — `fp-contract=fast` → NaN in the F distribution — reproduced, three defects fixed

The reporter's symptom is **real** and I reproduced it, but the
diagnosis is wrong and so is the work-around. He blames
`gsl_sf_beta_inc()` in `specfunc/beta_inc.c` and names
`gsl_sf_beta_inc_AXPY`, neither of which exists in `specfunc`. The
culprit is a *second*, independent copy of the same continued fraction in
`cdf/beta_inc.c`, reached through `beta_inc_AXPY()` which `cdf/fdist.c`
includes.

**Reproduction** — the machine has FMA, so the reporter's condition is
directly testable. Built the library twice with `gcc -mfma`, once
`-ffp-contract=fast` and once `-ffp-contract=off`:

```
=========== -ffp-contract=off ===========
0 NaN out of 1260 grid points
reported case: 0.051660470639181395
=========== -ffp-contract=fast ===========
0 NaN out of 1260 grid points
reported case: nan
```

`specfunc/beta_inc.c` is innocent — `gsl_sf_beta_inc_e` returns the
correct value in both builds. `cdf/beta_inc.c` has three faults:

**1. `max_iter = 512` is too small.** Measured by instrumenting the loop
and raising the limit:

```
  max_iter   iters used            cf          last |d-1|
       512          35    22098.9873275606         0     <- contract=off
       512         512    22098.9873275925  8.66e-15     <- contract=fast, exhausted
      1024         532    22098.9873275921  3.33e-16     <- converges at 532
```

Contraction rounds the terms differently and the convergence test stalls;
**35 iterations becomes 532**, twenty over the limit. `beta_cont_frac`
then returned `GSL_NAN`, and the caller had no error code to inspect
because `gsl_cdf_fdist_Q` returns a plain `double`.

**2. `GSL_NAN` used as a sentinel, then divided by.** This is the deeper
fault and the one that makes the NaN *survive*:

```c
if (fabs (den_term) < cutoff) den_term = GSL_NAN;
den_term = 1.0 / den_term;      /* NaN/NaN */
```

`specfunc/beta_inc.c` has the identical guard but clamps to `cutoff`. The
`cdf` copy poisons the term and then reciprocates it, so any trip of that
guard is an unconditional NaN regardless of iteration count. Latent
before; now clamped to match.

**3. `epsabs` cancels to exactly zero.** In `beta_inc_AXPY`:

```c
epsabs = fabs ((A + Y) / (A * prefactor / b)) * GSL_DBL_EPSILON;
```

For `fdist_Q`'s Q tail `A = -1`, `Y = +1`, so `A + Y` is **exactly zero**
and `epsabs = 0`. That disables the absolute convergence test
(`cf * |delta_frac - 1| < epsabs` can never fire), leaving only the
relative one — which needs many more iterations when `cf` is large, and
`cf = 22099` here. Confirmed by instrumenting: with `epsabs = 0` the loop
needs 532 iterations under contraction; with `epsabs` floored at one
epsilon the *relative* test settles it. Floored now.

I initially misdiagnosed this one, twice: I first assumed the values were
fine and only the error estimates were NaN (wrong — the values were NaN),
then I fixed `epsabs` alone and rebuilt, and it was *still* NaN, which is
what exposed the `max_iter` limit as the dominant factor. Three defects,
all needed.

**Verification.** After the fix, both builds give finite values agreeing
to 11 significant figures — `0.051660470639181395` vs
`0.051660470639254871` — which is exactly what a different rounding path
should produce. A sweep of 14400 points over the F, t and beta
distributions has no non-finite results in either build.

Negative control: rebuilding the contracted library with `max_iter` put
back to 512 makes the new test fail with

```
FAIL: gsl_cdf_fdist_Q(3.786820954867802, 1, 100000)
      (nan observed vs 0.0516604706391813953 expected)
```

so the test genuinely detects the bug rather than merely passing. Test
values for the three large-`nu2` cases are `1 - P` of the gp-pari values
already in `cdf/test.c`, so they are independent of the code under test;
GSL agrees with them to ~3e-15 against a `TEST_TOL6` of 2.3e-10.

Files: `cdf/beta_inc.c`, `cdf/test.c`. Commits 1b1d94ee9, 98e966458.

The reporter's workaround, building with `-ffp-contract=off`, is not
needed after this and should not be used as a substitute.

### `#57978` — cos_pi / sin_pi at infinity and NaN — applied

What the built library did before the fix:

    modf( inf) -> fracx=0 intx=inf
    modf( nan) -> fracx=nan intx=nan
    sin_pi_e( inf) -> status=0 val=0  err=0
    cos_pi_e( inf) -> status=0 val=1  err=0     <-- the bug
    sin_pi_e( nan) -> status=0 val=nan err=nan
    cos_pi_e( nan) -> status=0 val=nan err=nan

The reframing that settled the choice of return value: the `TWOBIG`
branches are *deliberate* large-argument approximations, not accidents.

    #define TWOBIG (2.0 / GSL_DBL_EPSILON)
    if(fabs(intx) >= TWOBIG) return GSL_SUCCESS;                 /* sin_pi */
    if(fabs(intx) >= TWOBIG) { result->val = 1.0; return GSL_SUCCESS; }  /* cos_pi */

Infinity is absorbed only because `fabs(inf) >= TWOBIG` holds, so this
is a non-finite argument falling into a finite-argument path — a domain
violation, not an accuracy problem.  Evidence for `GSL_EDOM` over the
patch's `GSL_ELOSS`:

* across `doc/*.rst`, 48 functions list `GSL_EDOM` alone, 3 list
  `GSL_ELOSS` alone (`polar_to_rect`, `angle_restrict_symm`,
  `angle_restrict_pos` — all pure accuracy cases) and 3 list both;
* `specfunc/error.h` fixes `DOMAIN_ERROR` as `NaN`/`NaN`/`GSL_EDOM`,
  and every domain-violation test asserts exactly that triple, e.g.
  `TEST_SF(s, gsl_sf_coupling_3j_e, (-1, 1, 2, 1, -1, 0, &r), GSL_NAN, GSL_NAN, GSL_EDOM)`;
* the patch's own `val = 0.0, err = 1.0` was not the value `cos_pi` was
  returning (1.0), so it was not a well-formed accuracy report either.

NaN was left propagating.  No specfunc test anywhere passes a non-finite
argument, and nothing treats NaN as a domain error — it flows through
the arithmetic and, where a convergence check exists, surfaces as
`GSL_EMAXITER` (`test_sf.c:600`, `gsl_sf_ellint_Kcomp_e(GSL_NAN, ...)`).

Blast radius: the only internal callers are `gsl_sf_bessel_Jnu_e` and
`gsl_sf_bessel_Yn_e`, reached only for `nu < 0.0`, propagating with
`GSL_ERROR_SELECT_4`.  The new guard can only fire when `nu` is itself
infinite, which was already meaningless.

Implementation notes: `sincos_pi.c` needed `#include <gsl/gsl_sys.h>`
(for `gsl_isinf`, which is what `config.h` maps `isinf` to) and
`#include "error.h"` (for `DOMAIN_ERROR`), following `psi.c`'s pattern.
Verified output after the change:

    sin_pi_e(   inf) -> EDOM  val=nan err=nan
    sin_pi_e(  -inf) -> EDOM  val=nan err=nan
    cos_pi_e(   inf) -> EDOM  val=nan err=nan
    cos_pi_e(  -inf) -> EDOM  val=nan err=nan
    sin_pi_e( 1e+16) -> SUCCESS val=0 err=0     (unchanged)
    cos_pi_e( 1e+16) -> SUCCESS val=1 err=0     (unchanged)

`sin_pi` and `cos_pi` had no `.. function::` entry at all, which is the
root reason this bug was ambiguous; both are now documented with
`.. Domain: -infinity < x < infinity` and
`.. Exceptional Return Values: GSL_EDOM`.  `gsl_sf_sin` and
`gsl_sf_cos` are also undocumented, with empty `Exceptional Return
Values:` lines — left alone as out of scope.


### `#57979` + `#39372` — non-finite arguments to the two `hypot` families — two independent defects fixed

The two reports were filed four months apart but name different
functions; they are grouped here only because they share a cause pattern.
Commit `05221394f` (2011) gave `sys/gsl_hypot` an infinity guard, but the
same treatment was never applied to the specfunc copy or to `hypot3`.

**`#57979` — `gsl_sf_hypot_e` (`specfunc/trig.c`).**  It selects
`min`/`max` with `GSL_MIN_DBL`/`GSL_MAX_DBL`, i.e. `(a) < (b) ? (a) : (b)`
and its `>` twin.  A NaN compares false against everything, so it is never
selected and both aliases land on the same argument.  Measured against the
built DLL with the default handler off:

    gsl_sf_hypot_e(NaN, 1.0) -> st 0  val 1.4142135623730951  err 6.28e-16
    gsl_sf_hypot_e(1.0, NaN) -> st 16 val inf                 err inf
    gsl_sf_hypot_e(inf, 1.0) -> st 16 val inf                 err inf

The first is the report: `min = max = 1.0` gives `sqrt(2)`.  The other two
reach the overflow test because `inf < DBL_MAX/root_term` is false for an
infinite `max`, so the routine reports `GSL_EOVRFLW` for a result the C99
specification says must be `+Inf` with no error.  Fix: an infinite argument
returns `+Inf` with `err = 0`, otherwise a NaN propagates as NaN, both at
`GSL_SUCCESS`.  This is the same shape as the applied `#57978`
(`sin_pi`/`cos_pi`): a non-finite argument falling into a finite-argument
path.  No blanket `isnan` guard is added — the two checks are the specific
exception handling the standard defines for `hypot`.

**`#39372` — `gsl_hypot3` (`sys/hypot.c`).**  It forms
`w = max(|x|,|y|,|z|)` and divides the magnitudes by it, so an infinite
argument gives `inf/inf = NaN`:

    gsl_hypot3(inf, 1, 1)   -> nan     (reported symptom)
    gsl_hypot3(inf, NaN, NaN) -> nan   (must be +Inf)

The top-level NaN propagation happened to work (`NaN/finite` is NaN), but
only accidentally.  The reporter is the GSL maintainer; his suggested
`gsl_hypot(gsl_hypot(x,y),z)` could not be taken as written, because
`gsl_hypot` raised the same overflow error that `#57979` is about.  Fix:
`gsl_hypot3` returns `+Inf` for any infinite argument, matching
`gsl_hypot`, and otherwise propagates NaN.

Independent reference: the C99 F.10.4.3 / POSIX behaviour of `hypot`
(`+Inf` if any argument is infinite even when another is NaN; otherwise
NaN propagates; no range error from finite inputs).  Verified by running
the built DLL.

Negative controls: reverting the `gsl_hypot3` guard makes `sys_test`
report five failures (all new); reverting the `gsl_sf_hypot_e` guard makes
`specfunc_test` report six (all new).  With both applied `sys_test` is
341/341, `specfunc_test` is clean and `ctest` is 56/56 on MSVC.

Documentation: the infinity rule was added to `gsl_hypot`/`gsl_hypot3`
(`doc/math.rst`, `doc_texinfo/math.texi`) and to
`gsl_sf_hypot`/`gsl_sf_hypot_e`, which previously had an empty
`.. Exceptional Return Values:` line, now `none`
(`doc/specfunc-trig.rst`, `doc_texinfo/specfunc-trig.texi`).

Commits c1df353ae (#39372) and 1a470222a (#57979).


### `#66808` (+ `#52359`, `#52570`, `#51000`, `#58067`) — Airy functions — three separate defects fixed

The umbrella report is "Ai is inaccurate away from the origin"; the
attached file adds integer values and derivative-at-zero vectors.  I
measured GSL against references computed from the Maclaurin series in
200-digit `decimal`: the *value* error in the oscillatory region is
`~eps*|x|^{3/2}` (2.2e-14 at `x = -10`, ~100 ulp) and 1-3e-15 for `x > 0`.
That is the phase, not the series: `gsl_sf_cos_e` is good to 1e-16 up to
`|t| ~ 1e5`, and the modulus is good to 1e-16.

Three real bugs sat underneath it, all in the `x < -1` branch:

* `#52359` — `mod.err`/`phase.err` divided by the Chebyshev corrections,
  which can be exactly zero (`x = -1.842761151977744`), giving `err = inf`.
  The formula was also wrong by a factor `2m/|rm|` elsewhere.  Rewritten
  from the full expressions.
* `#52570` / `#51000` — for `|x|` above ~1e13 `gsl_sf_cos_e` cannot reduce
  the phase and returned values outside `[-1,1]`, then `inf`, so
  `Ai(-1.14e34)` was `-inf` and `Ai'(-5.6e102)` had the wrong sign, with
  status `SUCCESS`.  Now: zero value, modulus/amplitude as the error,
  `GSL_ELOSS`.
* The derivative branch propagated `|val|*p.err` instead of `|ampl|*p.err`.

`#58067` (`Ai(113)` → `GSL_EUNDRFLW`) is not a bug: the value genuinely
underflows to zero; it is now documented.

**The accuracy itself was not improved, deliberately.**  A double-double
(Dekker) phase prototype (split `p = -0.625 + rp`, `cos(hi) - lo*sin(hi)`)
recovers the arithmetic rounding, but the error floors at
`eps*|rp|*|x|^{3/2}`, the precision of the double SLATEC coefficient
tables — a 5-10x gain at best, short of `TEST_TOL0`.  Recorded in
FORKNEWS; the added vectors use `TEST_TOL2` for `x < 0`.

**Fallout:** correcting the Airy error bars exposed an optimistic error in
`bessel_olver.c` (the Olver argument's rounding was never propagated);
fixed using the Airy ODE `Ai'' = arg*Ai`.  This is why three Bessel
vectors at extreme arguments are involved.  With the Airy change reverted
`specfunc_test` reports 55 failures, all in the new vectors.

Files: `specfunc/airy.c`, `specfunc/airy_der.c`, `specfunc/test_airy.c`,
`specfunc/bessel_olver.c`, `doc/specfunc-airy.rst`.
Commits `ac5f72d98`, `790fe1078`, `86cb80782`.


### `#42472` — `gsl_linalg_HH_solve` / `HH_svx` — three defects, fixed by a tall QR

The report has two parts, and the posted patch only addresses the first.

1. **`HH_svx` dimensions transposed.**  `hh.c` rejected `A->size1 > A->size2`
   as "underdetermined" (backwards - that is overdetermined) and set
   `N = size1`, `M = size2`, so the Householder loops ran over the wrong
   dimension.  For an `M x N` system with `M > N` the routine bailed out
   immediately.

2. **`HH_solve` overflow.**  It checked only `A->size2 == x->size` and then
   did `gsl_vector_memcpy (x, b)`; with `b` of length `M > N` this writes
   past `x`.  The reporter explicitly left this to the maintainer.

3. **The posted patch is necessary but not sufficient.**  I checked it on the
   report's own 5x3 example in the actual C (not by reading): swapping the two
   aliases makes the guard pass, but the body is still the *square*
   [Engeln-Mullges + Uhlig, Alg. 4.42] reduction, and it does not compute a
   least-squares solution for a tall matrix.  So the fix replaces the body
   with a proper tall Householder QR, shared by both entry points:
   `HH_solve` requires `b` of length `M`, `x` of length `N`, and returns the
   least-squares solution when `M > N`; `HH_svx` keeps the in-place contract
   and therefore requires a square matrix.

The existing tests only built `dim x dim` data, which is why the transposition
survived this long.  The new vectors include the report's polynomial fit
(reference from Burkardt's `QR_SOLVE`), rectangular Cauchy and Vandermonde
systems with exact solutions, and the error paths.  With the old `hh.c`
`linalg_test` reports 10 failures, all in the new vectors.

Files: `linalg/hh.c`, `linalg/test.c`, `doc/linalg.rst`,
`doc_texinfo/linalg.texi`.  Commit `1f84bed77`.


### `#43809` (+ `#28267`) — `gsl_sf_hyperg_1F1` for negative `a`, large `x` — fixed; `#28267` still open

Both reports carry only a test program as the attachment (`gsl_hyperg.c`,
`hyperg1F1.c`), not a proposed patch.

`gsl_sf_hyperg_1F1_e` routes `a < 0, b > 0` through the Kummer
transformation `M(a,b,x) = exp(x) M(b-a,b,-x)` and evaluates the
transformed call with `hyperg_1F1_ab_pos()`.  For `x > 0` that call has a
negative argument and `b < b-a`, so it takes the "recurse down in b"
branch.  Instrumenting it showed an accurate seed at `b ~ a` collapsing
to `1.2e-25` after 32 backward steps, where the true value is `2.72e-45`;
the recurrence is being run in the unstable direction for these
parameters.  The alternative `a0` branch is only valid for
`a > 0.5(b-x)`, and the large-`abs(x)` asymptotic only engages for
`x >~ abs(a)*abs(1+a-b) ~ 1118`, so the whole range `x` in roughly
`[66, 1600]` used the unstable branch.  For the reported parameters
(`a = 1 - 32.950611846591684`, `b = 2`, `x = 242.7876`) GSL returned
`3.39e80`; the true value is `7.51e60`.

The direct Taylor series in the original `(a,b,x)` is well conditioned
for large `x` here (no cancellation at the reported point), so the fork
evaluates it as well and uses it when its own error estimate is below
`1e-10` (at least ten significant digits).  A direct comparison of the
two error estimates is *not* usable: the Kummer estimate is wildly
pessimistic (correct values reported with `err/abs(val)` of order 10), so
"smaller estimate wins" gave 711 regressions over the test grid; the
series-only criterion gives none.

Independent reference: mpmath `hyp1f1` at 60 digits.  Over a 2400-point
grid (`a` in `[-50,-1]`, `b` in `[0.5,5]`, `x` in `[1,700]`) the fix
rescues 345 points and makes none worse.  The regression test uses the
bug's own parameters and fails on the pre-fix library with `3.39e80`;
full suite 56/56.  Files: `specfunc/hyperg_1F1.c`,
`specfunc/test_hyperg.c`.  Commit `52505315d`.

`#28267` is the same defect.  Its `(-26.1, 2, 100)` vector improves from
2.5% to `3.1e-10`, but `(-37.8, 2.01, 103.58)` stays at ~2%: in the
transition region `x ~ abs(a)^2` the series loses too many digits and its
estimate is too poor to be trusted, so the fix deliberately does not use
it there.  `#28267` remains open.


### `#39057` - `gsl_cdf_chisq_Pinv` "fails for some values" - fixed, but the report's expected value is wrong

The call is real: `gsl_cdf_chisq_Pinv (0.5, 0.01)` aborts with "inverse
failed to converge".  The report's expected value,
`0.99477710813146`, is not the inverse at all - it is
`gsl_cdf_chisq_P (0.5, 0.01)`, the forward CDF at `x = 0.5`.  The
inverse is `7.016667765235591e-61` (scipy 1.15 and mpmath at 60 digits
agree).  The test carried this wrong constant and had been disabled.

The cause is that the initial-approximation branching in
`gsl_cdf_gamma_Pinv` is only valid for shape `a >= 1`.  For the reported
`a = 0.005` it starts at `x0 = a` and the damped fixed-step iteration
halves `x` on each step; the true median is `3.51e-61`, about 200 steps
away, so it fails at the 32-step limit.  `gsl_cdf_gamma_Qinv` has the
same defect, and no convergence check, so it emitted wrong values
silently: for `a = 0.65` it returned `246.7` where the root is `500`.
The generated `test_auto` vector
`gsl_cdf_chisq_Qinv (5.8402405187288964e-219, 1.3)` expected `1000` and
got `493.4`, and now passes.

The solver is reworked as a bracketed Pegasus (regula falsi) iteration
on `t = log(x)`; the signed residual is monotone in `t` with limits
`-target` and `1 - target`, so a bracket always exists in the
representable range.  Both entry points delegate to whichever tail is
no larger than one half, and the forward incomplete gamma is evaluated
through the status-returning routines, with the complementary
continued fraction where the large-x expansion of `Q` would raise.

The fix was checked against `scipy.special.gammaincinv` over shapes
`1e-3 .. 1e6` and a wide range of tail probabilities.  Remaining
differences are only subnormal roots (unrepresentable, returned as `0`)
and points where GSL's forward `gamma_Q` is itself inaccurate (for
`a = 0.001`, `Q = 1e-12`, the forward value is off by about `3e-4`).
Negative control: with `cdf/gammainv.c` reverted, the new `cdf/test.c`
vectors abort at `gammainv.c:111`.  Full suite 56/56.  Files:
`cdf/gammainv.c`, `cdf/test.c`.  Commit `9befdae95`.


### `#53451` — elliptic `Pi` for a negative characteristic — the missing piece was the Cauchy principal value

The title names `gsl_sf_ellint_Pcomp(k, n, mode)` and says it returns
`NaN` for `mode < -1`, but the numbers in the report are the
*incomplete* `gsl_sf_ellint_P(phi, k, n)`:

    P_e(1.3, 0.5, -0.1) = 1.435170...    = EllipticPi(0.1, 1.3, 0.25)
    P_e(1.3, 0.5, -1.1) = NaN            vs EllipticPi(1.1, 1.3, 0.25) = 4.72112

and the convention matches once GSL's `n` is taken as the negative of
A&S's (as the manual states).  `Pcomp` has the same defect, just with
`p = 1 + n`.

The root cause is not a bug in the Pi formulas: both call
`RJ(x, y, 1, p)` with `p = 1 + n sin^2(phi)`, which is negative exactly
when the integrand has a pole in `(0, phi]`, i.e. `n < -csc^2(phi)`.
The Riemann integral diverges there, but the Cauchy principal value is
well defined and is what Mathematica, Boost.Math and mpmath's analytic
continuation (real part) return.  `gsl_sf_ellint_RJ_e` simply rejected
`p < 0`, so GSL returned `NaN`.

Applied: `gsl_sf_ellint_RJ_e` now uses the DLMF 19.20.14 transformation
(the one Boost.Math implements) for `p < 0`, so `P`, `Pcomp` and direct
`RJ` callers all get the principal value.  `p = 0` remains a domain
error.  A pre-existing error-estimate defect in `gsl_sf_ellint_P_e` was
fixed at the same time: it propagated the terms carrying the
characteristic with a signed `n/3`, which for `n < 0` could make the
reported error negative.

The formula was validated against `mpmath.elliprj` (40 digits) with the
same double arguments: 2.6e-15 worst relative difference over 25 random
`p < 0` cases.  Note the conditioning: the exact-real value differs
from the double-argument value by ~1e-14 near the pole, which is why
the regression vector for the reported point is quoted at `TEST_TOL2`.

Negative control: with the new RJ branch disabled, four of the five new
vectors fail; the fifth (`n = -0.1`) has no pole and guards the error
estimate.  Full `ctest` 56/56 on MSVC.  Files: `specfunc/ellint.c`,
`specfunc/test_sf.c`, `doc/specfunc-ellint.rst`,
`doc_texinfo/specfunc-ellint.texi`.  Commit `2673ce0d8`.


## Rejected

### `#39292` — missing `0.5` factor in `coulomb.c` — rejected, and the test gap is now closed

**Resolved by commit `192ec9cbf`**: ten test vectors were added to
`specfunc/test_coulomb.c` covering the Steed branch, so the patch can no
longer pass silently. The rest of this entry is the reasoning.

**Follow-up (commit `6eef17cf3`)**: two of those Steed vectors were still
held at `TEST_TOL2` even though their own comment records that the
independent reference is only good to `6e-14` relative there, and they
failed on arm64 macOS. Their tolerance is now `TEST_TOL3`; the half-`C`
patch is still detected.

Posted change at `specfunc/coulomb.c:1124`:

```c
- const double C = sqrt(1.0 + 4.0*x*(x-2.0*eta));
+ const double C = 0.5 * sqrt(1.0 + 4.0*x*(x-2.0*eta));
  const int N = ceil(lam_F - C + 0.5);
  const double lam_0   = lam_F - GSL_MAX(N, 0);
```

**This makes things worse, and I only found out by checking the geometry
rather than by reading the report.** My first instinct was that the
report was right about the quadratic — it is, up to a point — and that
the fix needed care. It needed more than care: the factor is wrong.

Write `Q = x(x - 2 eta) > 0`, which is guaranteed in this branch since
`x > 2 eta`. The turning point of the WKB region is where

    rho_ghalf = sqrt(x(2 eta - x) + lam(lam+1))

vanishes, i.e. `lam(lam+1) = Q`, so

    lam_turn = 0.5*(sqrt(1 + 4Q) - 1)        hence  sqrt(1 + 4Q) = 2 lam_turn + 1

`coulomb_jwkb()` needs `lam > lam_turn`. The code picks `lam_0` by
rounding *down* in whole units of lambda, because the recurrence
advances one lambda at a time. Rounding down needs the margin to be a
whole number of lambdas, which is exactly what the current expression
supplies:

| expression | value | margin above `lam_turn` |
|---|---|---|
| current `sqrt(1+4Q)` | `2 lam_turn + 1` | **1 lambda — correct** |
| posted `0.5*sqrt(1+4Q)` | `lam_turn + 0.5` | half a lambda — too small |

With a half-lambda margin, `ceil(lam_F - C + 0.5)` rounds *up* past the
turning point and `lam_0` lands **below** it. Measured over a grid of
`(eta, x)` with `lam_F` sampled above the turning point:

```
       Q  lam_turn | current C = sqrt(1+4Q)     | patched C = 0.5*sqrt(1+4Q)
    0.25    0.2071 | lam_0= 0.5771 ok           | lam_0=-0.4229 BELOW TURN
       1    0.6180 | lam_0= 0.9880 ok           | lam_0=-0.0120 BELOW TURN
       6    2.0000 | lam_0= 3.3700 ok           | lam_0= 1.3700 BELOW TURN
      20    4.0000 | lam_0= 5.3700 ok           | lam_0= 3.3700 BELOW TURN
      99    9.4624 | lam_0=10.8324 ok           | lam_0= 8.8324 BELOW TURN
     400   19.5062 | lam_0=20.8762 ok           | lam_0=18.8762 BELOW TURN
    1e+04   99.5012 | lam_0=100.8712 ok          | lam_0=98.8712 BELOW TURN
```

Below the turning point `rho_ghalf` is imaginary and `coulomb_jwkb()`
cannot be used at all, so the patch breaks the branch it intends to
repair. The reporter's own claim — "none of the current test cases are
influenced" — is true, and I confirmed it: applying the patch leaves
`specfunc_test` at zero failures. That is the danger. The existing
vectors do not reach this branch with `lam_F` close enough to
`lam_turn` for the difference to show, so **the test suite cannot catch
this patch, and applying it would have silently degraded
`gsl_sf_coulomb_wave_FG_e`**.

Independent confirmation that the current code is right: in the
`x > 2 eta` branch the library satisfies the Coulomb Wronskian
`F Gp - G Fp = -1` to the last bit across the whole sampled range,
which it could not do if `lam_0` were being placed wrongly.

```
  eta      x    lam   F*Gp - G*Fp    expected     rel dev
    0    2.5      0              -1            1            0
    0      10      0              -1            1            0
  0.5      20      0              -1            1            0
   -2      50      0              -1            1            0
```

Not applied. Note the submitter of the original report is Alexey A.
Illarionov, who also wrote the 3j/6j patch accepted as #39473, so this
is not a careless report — but the `+0.5` in `ceil(lam_F - C + 0.5)` is
precisely what makes the full `sqrt` correct, and halving `C` breaks the
margin the rounding depends on.

**How wrong, measured.** Applying the patch and sweeping the branch,
15 of 789 sampled points move by more than 1e-9 in `F`:

```
   eta    lam        x           dF          dFp           dG          dGp
   -1     30       12     1.54e-06     1.54e-06     1.54e-06     1.54e-06
   -2     30       20            1            1     2.13e+05     4.68e+05
   -2     30       35       1.01            1          175          758
   -2     30       60       0.601        5.31           8.4        0.699
   -2     30      100       0.033        4.36          4.75       0.0452
```

`F` is wrong by a factor of two, or becomes entirely spurious, and `G` is
off by five orders of magnitude, at points the suite never visited. This
is a much stronger statement than my first pass at it, which rested on
`lam_0 < lam_turn` reasoning that turned out to describe a *different*
regime from the one the suite reaches.

**The test gap, and what I got wrong filling it.** My first attempt at
reference values was wrong in three successive ways, each caught by a check
rather than by reading:

1. Integrating the Coulomb ODE from `x0 = 1e-6` in `x`. The `1/x²` term is
   stiff there; the initial slope came out as 0.54 instead of 1 because I
   had written `sp += k b_k (L+1) pw/x0` where the exponent of `pw` was
   off by one. Fixed by evaluating `powl` per term.
2. Then the RK4 converged at second order, which I misread as "needs more
   steps" and duly added — making it *worse*. The stage-3 state was being
   advanced by stage-1 derivatives' half-step instead of its own. Once
   fixed, the error ratios matched `(N_{k+1}/N_k)^4` and it was simply a
   matter of step count.
3. The reference then agreed with the library to 1e-15 at `eta = 0` but
   disagreed by 2–6% at `eta != 0`, and I initially read that as the
   library being wrong. It was my `Clam`: GSL's `C_lam` uses the
   *magnitude* of `Gamma(lam+1+i eta)`, so the logarithm is the real part
   of `lngamma` of the complex argument, not `lgamma(lam+1)` on the real
   axis. With that corrected the reference matched to 1e-15 everywhere.

Two derivations I had to discard as wrong rather than merely unhelpful: the
Bessel order is `lam + 1/2`, not `sqrt(lam(lam+1) + 1/4)`; and the Wronskian
of the *reduced* pair is proportional to `1/x`, so it pins `B_lam` only
in combination with the small-`x` condition, not alone. Measuring the
constants directly settled both: `A_lam = B_lam = ±sqrt(pi/2)`,
independent of `lam` to 1e-15.

**A negative result worth keeping.** The Wronskian
`F G' - G F' = -1` is exact to 1e-16 in *both* builds, so it cannot detect
this patch. It is a good check on the values but useless as a regression
test here. Only comparing against an independent reference detects it.

**And a tolerance lesson.** `test_sf_check_result()` has a second
condition beyond the relative tolerance:

```c
if (d > 2.0 * TEST_SIGMA * r.err) s |= TEST_SF_INCONS;
```

so an expected value must lie inside the library's own reported error bar.
At two of the four negative-`eta` points the library's error is 6e-14
relative while the reference is only good to 5e-13, so quoting the
reference there produces `TEST_SF_INCONS` — which would assert that the
library's error estimate is too small, a different claim from the one
these vectors are meant to make. The values are therefore quoted at the
library's full precision, with the independent agreement (1e-14 or better
at three of four points) recorded in the comment as the justification.
That is circular in form and not in substance, but it is worth being
explicit about, because the obvious alternative is a test that fails for
a reason unrelated to what it is testing.

**Follow-up (2026-10-02, commits `491f38efa`/`a71516e66`): the rejection
above is not supported, and the "factor of two" measurement is wrong.**
This surfaced while fixing the arm64 CI, which still failed at
`TEST_SF_INCONS` after `6eef17cf3`.

* Re-applying `0.5 * C` to the current tree moves the four test values by
  at most a few times `1e-13` (fracdiff `5.0e-13` at `eta = -2,
  lam_F = 30, x = 20`); it does **not** change `F` by a factor of two
  there. The values only move appreciably near the turning point
  (`7.5e-6` at `x = 15`, `3e-5` at `x = 12`).

* An independent integration of
  `y'' = [lam(lam+1)/x^2 + 2 eta/x - 1] y`, started from the Frobenius
  series `b_k = (2 eta b_{k-1} - b_{k-2}) / (k (2 lam + k + 1))` and
  normalised by `C_lam`, was validated against the exact `eta = 0`
  relation `F_l(0,x) = sqrt(pi/2) sqrt(x) J_{l+1/2}(x)` to `4e-14` at
  `l = 30, x = 12, 20, 35`. It agrees with the **patched** build:
  relative error at `eta = -2, lam_F = 30` is `3e-5` (unpatched) vs
  `2e-14` (patched) at `x = 12`, and `7.5e-6` vs `2e-14` at `x = 15`.
  A Kummer-series evaluation agrees at `x = 12`. So `0.5 * C` appears to
  *fix* an accuracy defect below the turning point, not introduce one.

* The Wronskian negative result above still holds: it cannot see the
  patch. The over-tight `TEST_TOL2` was what made the old negative control
  fire, not a real factor-of-two error.

**DONE (commit 90d9037a5).** `coulomb.c` now uses
`C = 0.5 * sqrt(1 + 4Q)` (equivalent to `lam_turn + 1/2`), a fifth Steed
vector at `eta = -2, lam_F = 30, x = 12` guards it (fracdiff `1.49e-5`
when the `0.5` is reverted), and the `#39292` verdict is reversed. The
arm64 CI fix (commits 491f38efa/a71516e66) still uses
`test_sf_check_val` with `TEST_TOL5` at the two cross-target-sensitive
points.

### `#64777` — complex LU "returns incorrect results" — not a bug, and demonstrably so

The reporter's own matrix is singular:

    A = [1 2 3]
        [4 5 6]
        [7 8 9]

    det(A)        = 6.66e-16          (i.e. zero)
    U diagonal    = 7.0000 0.8571 0.0000   <- last pivot is zero

The Julia output he quoted as "the correct LU decomposition" ends in a row
of zeros, which is the signature of a rank-deficient matrix, not a bug.
The two matrices he pasted are not an LU pair at all: they are the
pivoted matrix `PA` and the upper factor `U`, side by side, and they do
not multiply back to `A`. And his code allocates
`gsl_permutation_calloc(3)` and never reads it back.

Running his exact program and splitting the packed factor properly:

    real    LU permutation: 2 0 1   signum 1
    complex LU permutation: 2 0 1   signum 1     <- identical

    max |L*U - P*A| = 0.000e+00      the identity holds exactly

So `gsl_linalg_complex_LU_decomp` is correct, the real and complex paths
agree on real-valued input (they share the BLAS `izamax`/`idamax`
convention), and the only thing missing from his comparison was the
permutation. This matches Christian Krueger's reply on the tracker: LU
with partial pivoting is not unique, `PA = LU` holds with whichever
permutation GSL chose, and comparing raw factors from tools that pivot
differently is meaningless.

No code change. Recorded so the next person to read the report does not
re-investigate it.

### `#66862` — duplicated complex sin/cos "produce differing results" — not a bug

**The submitter withdrew it himself.** The follow-up comment says "the
findings were due to rounding errors and both implementations ultimately
rely on the identities" and, of the strongest claim, "**finding 3 can be
ignored**".

Nothing was left to fix. `gsl_complex_sin` (`complex/math.c:450`) and
`gsl_sf_complex_sin_e` (`specfunc/trig.c:362`) compute the same closed
form; they differ only in how `cosh`/`sinh` are obtained:

| | `\|z\| < 1` | `\|z\| >= 1` |
|---|---|---|
| `complex/math.c` | libm `cosh`, `sinh` | libm |
| `specfunc/trig.c` | `sinh_series`, `cosh_m1_series` | `0.5*(exp(z) +- 1/exp(z))` |

The specfunc series path is the *more* accurate branch. Measured worst
relative error against a long-double reference, over the whole grid:

| `zi` | sin | cos | status |
|---|---|---|---|
| 0 | 4.5e-17 | 8.8e-17 | SUCCESS |
| ±0.5 | 1.4e-16 | 9.9e-17 | SUCCESS |
| ±0.9 | 1.8e-16 | 1.5e-16 | SUCCESS |
| ±1 | 1.6e-16 | 1.3e-16 | SUCCESS |
| ±5 | 1.2e-16 | 1.8e-16 | SUCCESS |
| ±30 | 2.3e-16 | 2.3e-16 | SUCCESS |
| ±100 | 2.3e-16 | 2.3e-16 | SUCCESS |
| ±700 | 1.2e-16 | 1.2e-16 | SUCCESS |
| ±710 | inf | inf | `GSL_EOVRFLW` (16) |

Every value is at or near machine epsilon, and the overflow past
`GSL_LOG_DBL_MAX` is reported rather than returned silently.

**Finding 3 was arithmetically void.** He claimed
`cos(-6.2-5.3i)` should equal `cos(6.2-5.3i)` because cosine is even.
But negating `-6.2-5.3i` gives `6.2+5.3i`, so that pair is not
`cos(z)`/`cos(-z)` at all. Evaluating the identity at his four points:

    -6.2  -5.3i ->  99.824520   8.322726i
     6.2  -5.3i ->  99.824520  -8.322726i
    -6.2   5.3i ->  99.824520  -8.322726i
     6.2   5.3i ->  99.824520   8.322726i

consistent — the imaginary part flips with `I` alone, which is correct.
The values he quoted (`98.008630 + 9.833157i`) match neither
implementation nor the identity, so they came from an untrustworthy
intermediate.

**The posted test patch was unusable as-is**, and this is worth
remembering as a pattern: of the 22 `TEST_SF_2` lines in
`specfunc_test_sf.diff`, **3 were exact duplicates** (`sin_e(1.0, 0.0)`,
`sin_e(1.0, 1e-9)`, `sin_e(1.0, 5.0)` each appear twice — evidently two
pasted runs) and **13 were commented out** with the note `// accuracy`.
Every surviving line used `TEST_TOL0` on values up to 169.7, demanding
`2*eps` relative accuracy for results built on `exp` of a large
argument. His `test_complex.c` is a standalone program with its own
`main`, not wired into any build.

What was taken from it: the coverage gap, rewritten. Recorded in
FORKNEWS as `[rejected]` (commit 758c15422) and the vectors added in
48a49003f.

Also noted and **not** done, as new API: no non-`_e`
`gsl_sf_complex_sin`/`gsl_sf_complex_cos` exists, no
`gsl_sf_complex_logcos_e` to match `gsl_sf_complex_logsin_e`, and
`lnsin` is inconsistent with `gsl_sf_lngamma`/`gsl_sf_lnbeta`. His
suggested delegation refactor was rejected: it would change results for
a module with no error-estimate machinery, and with the discrepancy
attributed to rounding there is no correctness case for it.

### `#52321` — bidiag_unpack2 householder call — the patch is wrong

Posted change: `gsl_linalg_householder_hm` → `gsl_linalg_householder_mh`
in `gsl_linalg_bidiag_unpack2()`.  Applied cleanly, but `unpack2` is not
broken.  Measured against `gsl_linalg_bidiag_unpack()` on random 6x4
matrices:

    max |V_unpack - V_unpack2|      = 0.000e+00
    max |U_unpack - U_unpack2|      = 0.000e+00
    max |A - U B V^T| via unpack    = 1.776e-15
    max |A - U B V^T| via unpack2   = 1.776e-15

Both produce bit-identical `U` and `V`, and both satisfy `A = U B V^T`.
`unpack2` differs from `unpack` only in that it accumulates `U` in place
with `gsl_linalg_householder_hm1`, which is the correct routine for that
case: `hm1` treats the first column of its argument as the Householder
vector, and in the in-place accumulation that column is exactly the
compressed vector that `bidiag_decomp()` stored.  Applying the patch would
break working code.

The reasoning is easy to get wrong here because `unpack` and `unpack2`
look interchangeable in the source: `unpack` uses
`gsl_linalg_householder_left` where `unpack2` uses `..._hm`, and both
multiply from the left.  They agree.

### `#58066` — digamma domain test — **I was wrong; the patch was right, and it is now applied**

**Corrected and applied as commit `d53d00109`.** My first reading of this
was wrong and the reasoning below is kept only to show what the mistake was.

My original note read:

> Posted change replaces `x == 0.0 || x == -1.0 || x == -2.0` with
> `x <= 0.0 && round(x) == x`, which turns *every* negative integer into a
> `DOMAIN_ERROR`. The special case exists because `psi` has poles at the
> non-positive integers; broadening it to the whole set is not a fix.

The last sentence contradicts itself. `psi` has poles at *all* the
non-positive integers, so widening the check to exactly that set is the
correct fix, not a regression. My error was to read
`x == 0.0 || x == -1.0 || x == -2.0` as a deliberate list of the poles,
when in fact it is an incomplete one: the remaining integers were supposed
to be caught by the guard

    if(fabs(sin(M_PI*x)) < 2.0*GSL_SQRT_DBL_MIN) DOMAIN_ERROR(result);

which assumes `sin(M_PI*x)` underflows at an integer. It does not.
`sin(M_PI*(-3))` rounds to about -3.7e-16, fifteen orders of magnitude
above the `2*GSL_SQRT_DBL_MIN` of roughly 3e-154 that it is compared
against, so the guard never fires.

Measured, before the fix:

```
         x status                          psi
        -2 EDOM                            nan
        -3 SUCCESS       -8.55101692933585e+15   <- pole, should be EDOM
        -4 SUCCESS       -6.41326269700188e+15
      -100 SUCCESS        1.59927402055892e+15
```

Two more functions have the same defect, which is why the fix is one shared
predicate rather than a one-line change to `psi_x`:

* `gsl_sf_psi_1_e` (trigamma) catches -3 and -4 through the `fx == 0` test
  in its `x > -5` branch, but -5 and below go through the asymptotic branch
  and return about 2.6e31;
* `gsl_sf_complex_psi_e` rejects only the origin; `z = -1, -2, ...` return a
  large finite value because the reflection branch's `GSL_IS_REAL()` is just
  `gsl_finite()`, and `cot(pi z)` rounds to a large finite number instead of
  overflowing.

`gsl_sf_psi_n_e` for `n >= 2` is already correct, since it requires `x > 0`.

I used `floor(x) == x` rather than the posted `round(x) == x`; the two agree
here and `round()` is not used anywhere else in the library's live code.

The finding came from checking the domain while writing the
`gsl_sf_complex_psi_e` documentation (`#67689`), not from the tracker.
Non-integer arguments are unaffected: -1.5, -10.5, -100.5, -262144.5 and
the complex points still return their previous values. With the change
reverted `specfunc_test` reports 11 failures, all new vectors.

### `#59834` — movstat accumulator alignment — applied

**DONE (commit 4e4a88242).**  `mmacc_init()` places `state->maxque` at
`state->minque + deque_size(n+1)`.  A `deque` holds a pointer, so it
needs 8-byte alignment on the 64-bit targets; `deque_size(n+1) =
sizeof(deque) + 4*(n+1)` is 4 modulo 8 when `n` is even, so `maxque` is
misaligned.  Reproduced with clang `-fsanitize=alignment` (MSVC has no
UBSan and MinGW lacks `libubsan`, but LLVM clang 18 on this machine
works): `deque.c:58-61` reports "member access within misaligned address
... for type 'deque', which requires 8 byte alignment", matching the
report.

Auditing every accumulator for the same pattern also found `qnacc.c`:
`state->rbuf` follows a `5*n` int array, so its offset is
`sizeof(state) + 52*n`, which is 4 modulo 8 for odd `n`.  The sanitizer
reports it at `ringbuf.c:61`.  Both are fixed.

The earlier note here called the posted patch (#57859) a wrong fix; that
verdict was itself wrong.  Applied to a copy, the patch does align
`maxque` and its size accounting matches - it simply rounds the offset
up to a multiple of `sizeof(deque)`, and `sizeof(deque)` is always a
multiple of the alignment, so that is sufficient.  The note's stated
reason conflated alignment with `sizeof`, and `ringbuf_size(n)` is
`sizeof(ringbuf) + n*sizeof(double)`, already a multiple of 8, not
`n*sizeof(double)`.  The patch is still obscure and silently relies on
`minque` already being aligned; the fork instead adds a
`ringbuf_align()` helper in `ringbuf.c` and rounds each object offset
explicitly.

With the change reverted the sanitizer reports both misalignments; with
it applied the drivers are silent and the MSVC `ctest` suite is 56/56.

### `#47646` — `gsl_ran_beta` NaN for small parameters — code already upstream, test added

**DONE, test only (commit d64cc4d93).**  The reported defect is real but the
code fix landed upstream long ago as `05c5b5179` ("bug fix for #47646
(gsl_ran_beta)", 2016-09-15): below `a = 1` a gamma variate is exactly zero
most of the time (`gsl_ran_gamma(1e-5, 1)` gives exact zero for 198470 of
200000 draws), so the plain `x1 / (x1 + x2)` was `0/0` and returned NaN for
roughly 98 % of samples.  Upstream added a `(a <= 1 && b <= 1)` rejection
branch; this fork's `randist/beta.c` is identical to `savannah/master`.

The bug is still open upstream because no test came with the fix.  The
attachment `test_beta_small.c` (file 37915) works around this by computing
bin probabilities from `gsl_cdf_beta_P` via a new `testPDF_c`, but that
duplicates ~100 lines of `testPDF`, uses GSL's own CDF as the reference,
and carries debug `printf`s.  The other attachment, `beta_distribution.diff`
(file 37884), is a *reverse* of the upstream fix (it deletes the small
branch) and must not be applied.

Instead the fork adds `test_beta_small` to `randist/test.c` and drives it
through the existing `testMoments` harness.  With `a = b = 1e-5` the
distribution is symmetric about `1/2`, so `P(X < 1/2) = P(X > 1/2) = 1/2`;
the `a = b` symmetry is the reference and does not consult GSL's CDF.  The
intervals are `(-0.5, 0.5)` and `(0.5, 1.5)`, straddling the atoms at 0 and
1 that the strict inequalities in `testMoments` would otherwise discard
(a `(0, 1/2)` interval would see almost nothing).  The vectors are
registered last in `main` so the shared RNG stream is untouched for every
existing test - placing them earlier perturbed `dirichlet`'s mean test into
a failure, purely by shifting its variates.

Negative control: with the `(a <= 1 && b <= 1)` branch reverted to the
gamma ratio the two vectors fail (0.00759 and 0.00783 observed against 0.5
expected); with the fix they pass and the suite is 152/152, `ctest` 56/56.

File: `randist/test.c`.  See Savannah bug #47646.

### `#65912` — GSL_SET_COMPLEX — both posted variants are wrong

**DONE (commit 631da98f8).**  The native C11 branch built the value as
`(*(zp) = (x) + I*(y))`, a complex multiply-and-add, so signed zeros were
rounded away (`GSL_SET_COMPLEX(&z,-0.0,0.0)` gave `+0.0`) and a
non-finite part was corrupted (`(0.0, INFINITY)` gave `NaN`, from the
`0*inf` in `I*y`).  The macro now assigns the two components directly
through the same `GSL_REAL`/`GSL_IMAG` lvalues that `GSL_SET_REAL` and
the legacy fallback already use.  A block of exact checks was added to
`complex/test_source.c` (four fail under the native branch with the old
macro); `complex.test_c11` is now built by CMake so the native branch is
covered in CI (fork commit db5fcdf6d).  The rest of this entry is the
reasoning for the earlier rejection.

`GSL_SET_COMPLEX` is used with `float`, `double` and `long double`
pointers, and for the non-Annex-G fallback at
`gsl/gsl_complex.h:140`.  The `CMPLX(x, y)` variant silently truncates to
`double`; the `_Generic` variant covers the three complex types but drops
the `const` case and every non-complex fallback, so it will not compile
on MSVC or where `<complex.h>` lacks Annex G.  Would have to be written
from scratch.  Note the fork builds on MSVC, so this one matters here.

Confirmed while implementing: MinGW's `<complex.h>` defines `_Complex_I`
and `I` but **no** `CMPLX`/`CMPLXF`/`CMPLXL`, so both posted patches fail
to compile on a C11 toolchain that lacks Annex G `CMPLX`.

### `#61342` — test_c11 on ppc64 and sparc — not a library bug

No patch is attached.  The FAIL lines the reporter pasted come in adjacent
pairs with `real part = 0` and `imag part = 1`, which is the signature of
`gsl_complex_pow_real` being evaluated under a compiler mode that
mishandles the complex result.  The reporter's own follow-up points at the
`-ffast-math`-style options Gentoo and Fedora apply.  GSL does not support
compiling with `-ffast-math` (see `#39171`).  Nothing to change in the
library.

### `#50343` — mathieu_ce differs from Mathematica — not a bug

The submitter is the maintainer, and the thread settles it: Wolfram's
`MathieuCharacteristicA[0, -1]` picks a branch that Wolfram's own
documentation describes as arbitrary ("there is no general agreement on
how to define the branch cuts ... the implementation simply picks a
convenient sheet").  SPECFUN (via SciPy) and Scilab both agree with GSL's
`0.9975194` for `q = -1`, `z = 2 pi/180`.  No patch, no change.

### `#44865` — bsimp/msbdf `e5_bigt` is FMA-sensitive — patch not applied

The test `e5_bigt` in `ode-initval2/test.c` compares `bsimp` with `msbdf`
over `[0, 1e11]` using a pure relative tolerance of `1e-3`; with
`-ffp-contract=fast` the reported values are off by enough to fail.  The
posted change shortens the integration interval from `1e11` to `1e9`, which
makes the test pass by testing less, and that is not an acceptable fix.

Re-investigated 2026-10-02: the failure does not reproduce on the fork's
build hosts and the patch stays rejected.  At `t = 1e11` every component has
decayed to round-off (`y0 ~ 1e-49`, `y1, y2 ~ 1e-20`, `y3 ~ 1e-65`, from
`y0(0) = 1.76e-3`), so the relative comparison is between two independent
accumulations of rounding.  With FMA verified active - gcc 15 on x86-64 with
`-mfma -ffp-contract=fast` (`vfmadd` is emitted in `rhs_e5`) and MSVC with
`/fp:fast` - the two steppers agree to about `1e-4` relative and all 1151
tests pass.  The reported 18-22 % disagreement needs the submitter's
clang/PPC64 toolchain.  A test-quality change (an absolute-aware comparison)
would make the test robust, but it cannot be validated here and is not
required on any platform the fork builds and tests on, so it is left out.

### `#54925` — source_gemm_r loop reordering — out of scope

The patch swaps the `i`/`k` loop nesting in two of the four branches of
`cblas/source_gemm_r.h` to reduce cache misses.  It applies cleanly and is
plausibly a speed-up, but it is a performance change, not a bug fix or a
documentation change, and reordering the accumulation changes the floating
point result.  Outside the eligibility rule for this review.

**Superseded (commit `39cde03e2`).**  The "changes the floating point
result" claim is wrong: the patch only changes the traversal order, and
every `C(i,j)` still accumulates its terms in the same `k` order, so the
result is bit-for-bit identical.  Applied after the filter was widened;
see Group V at the end of this file.


### `#58763` — `gsl_root_fsolver_brent` "wrong results under valgrind" — not reproducible

Not a defect in the current source.  The reported output (first iteration
`root = 5`, bracket `[5, 5]`) follows from exactly one state: the Brent
third point `c`/`fc` being zero.  The same-sign fix-up at the top of
`brent_iterate` is then skipped, and the `|fc| < |fb|` swap collapses
`b = c = 5`, so `m = 0`, `fb = 0` and the routine returns immediately.
Brent is the only bracketing solver whose state carries a third point,
which is why only brent was affected.  `brent_init` initialises all eight
fields (`state->c = x_upper`, `state->fc = f_upper`), and `roots/brent.c`
is identical to the 2.5 and 2.6 releases.  A 2025 tracker comment also
failed to reproduce it with 2.6.

Verified on the fork 2026-10-02: the reporter's `demo.c`/`demo_fn.c`,
linked against the current shared library, prints the report's expected
(non-Valgrind) table - first iteration `[1.0000000, 5.0000000]`
`root 1.0000000`, converged on iteration 6 at `2.2360634`; `ctest -R roots`
passes.  The library the reporter linked must have had a `brent_init` that
omitted `c`/`fc`.  No change made.


## Deferred, pending a decision

### `#64549` — interpolation test cases

**Applied 2026-10-04** as one clean test rather than the four posted:
`test_eval_wrappers()` in `interpolation/test.c` (commit 7b6979126).
See `Applied` above and the `FORKNEWS` entry.

The posted patch adds four `test_hstaxe_1..4` functions to
`interpolation/test.c`.  Two of the four are near-duplicates of the
others, all four leak the `gsl_interp` and `gsl_interp_accel` they
allocate, two declare unused variables (`test_table`, `min_size`,
`wav_min`, `wav_max`), and `test_hstaxe_2` asserts the result is 0.0.
Low value as posted; the one clean test covers the plain wrappers and a
partial interval instead.

### `#64851` — gsl-config exit status

Superseded.  The one-line change (`usage; exit 1` → `usage 1`) is already
equivalent to what this fork's own `gsl-config` work does with
`usage 1` plus `exit $1`.


## Not yet reviewed

Everything on the curated list has now been reviewed; what remains are
feature-shaped patches excluded by choice (see the eligibility rule in the
discussion): `#66573`, `#66574`, `#66575` (complex SVD suite), `#66695`
(Feagin high-order ODE solvers, +5866/-35), `#66576`, `#67359`, `#66922`,
`#66886`, `#66850`, `#66842`, `#66834`, `#66742`, `#66767`,
`#66775`, `#66949`, `#66816`, `#68367`, `#60457`.

Everything else from the earlier list is in `Applied` or `Rejected` above.

The `source` rows (whole replacement files rather than diffs) are mostly
the 2025 specfunc/randist series from one contributor: `#66816`,
`#66834`, `#66842`, `#66850`, `#66886`,
`#66922`, `#66949`, `#67359`, `#67728`, `#66767`, `#66775`, `#66800`,
`#67621`, `#42472`.  These need the patch to be reconstructed by
hand; `#47345`, `#60371`, `#67621`, `#42472`, `#47646`, `#66862`, `#66877`
and `#66880` came out of this group and are now applied.


## Offline source: working without Savannah (2026-10-03)

`savannah.gnu.org` became unreachable during the review (the sweep was
rate-limited / the host refused connections), so the remaining bugs are
worked from copies held outside Savannah.  None of this is part of the
fork's record; the material lives under the git-ignored `temp/`.

Three sources were evaluated:

* **bug-gsl mbox archive** — `lists.gnu.org/archive/mbox/bug-gsl/`, one file
  per month, 2002-10 … 2026-09, 277 files, 82 MB.  Savannah mirrors every
  bug notification there as `[bug #NNNNN]`, so a bug's whole discussion is
  recoverable.  Downloaded in bulk (one request per month, small delay) to
  `temp/bug-gsl-mbox/`.  This is the primary fallback: it covers the bug
  *traffic*, which is what a review needs, even when the Savannah page is
  gone.
* **CNMAT `BUGS`** — `github.com/CNMAT/gsl/blob/master/BUGS`, a dump of the
  bug database generated 2011-04-01, whole bug text for everything up to that
  date.  Kept as `temp/cnmat_BUGS.txt`.
* **Wayback Machine** — the CDX listing
  `web.archive.org/cdx/search/cdx?url=savannah.gnu.org/bugs/&matchType=prefix`
  (one request) gives the snapshot id for every archived bug page.

Coverage of the remaining bugs is now complete: the subject parser was
fixed to accept a list tag before the bug tag (`[Bug-gsl] [bug #NNNNN]`),
which raised the captured `bug-gsl` threads from 119 to 215 of the 219 open
bugs, and every remaining bug has a merged dossier — the five earlier thought
to lack text (`50712`, `51104`, `52351`, `53903`, `53904`) turned out to be
present all along.  Text is merged into per-bug dossiers under
`temp/savannah-store/` (see its `README.md`), which is what the review now
reads instead of the tracker page.

The mbox host is not Savannah and is not rate-limited; keep the pulls to one
request per month.  `scripts/savannah_bugs.py` still needs Savannah for the
live list, so the 2026-10-01 cached `cache.json` is the metadata source until
it is reachable again; re-run `--refresh` once and diff the inventory when it
is.

A separate caution that this episode made concrete: the mbox archive also
proves that `git apply --check` on an attachment is not the only risk — the
*report itself* can be wrong, and on macOS `savannah.gnu.org` and the
`lists.gnu.org` archive are the two independent records that let a later
reader tell a withdrawn report from a real one (`#66862` is the example).


## Method notes

* **`-ffp-contract=fast` can be tested directly here.** This machine has
  FMA (`gcc -Q --help=target` shows `-mfma [disabled]`, i.e. available), so
  the `#64613` condition is reproducible without emulation. The reliable
  recipe, since `cmake -S . -B <dir>` + `ninja` works but linking against
  the DLL fails on the non-ASCII temp path:

      cmake -S . -B D:\gsl-fma-test\fast -G Ninja -DCMAKE_C_COMPILER=gcc `
            -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="-O2 -mfma -ffp-contract=fast"
      ninja -C D:\gsl-fma-test\fast gsl          # the gsl-randist/gsl-histogram
                                                # tools fail to link, but the
                                                # library is fine
      # archive the objects, because the import lib does not export everything
      Get-ChildItem ...\CMakeFiles\gsl.dir -Recurse -Filter *.obj |
        ForEach-Object { $_.FullName.Replace('\','/') } |
        Set-Content objs.txt -Encoding ascii
      ar rcs libgsl.a @objs.txt

  Three traps here, all of which look like a broken library: the response
  file needs **forward slashes** (backslashes get eaten and the linker
  reports `cannot find C:Users...`), `ar`/`gcc` on a response file need
  `-Encoding ascii`, and the DLL's import library does not list every
  symbol, so linking the test program against `libgsl.dll.a` gives
  `undefined reference to gsl_cdf_fdist_Q` even though the symbol is in
  the DLL. Archiving the objects sidesteps all three.

* Build: `build-cmake` is configured for MSVC, so a rebuild needs the
  developer environment:

      cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && cd /d D:\projects\crystal\gsl\build-cmake && ninja'

  gcc from MinGW is also installed but is not what that directory is
  configured for.  Baseline before any change: 56/56 `ctest` pass.

* Two of my own test programs gave wrong answers before they gave right
  ones, and in both cases the failure looked like a library bug:
  a `U2` matrix that was passed to a residual check but never populated
  from `A2`, and a `printf` that read `r.val` in the same argument list as
  the call that sets it (argument evaluation order is unspecified).  When
  a Savannah patch says the library is broken, check the harness first.


### `#36152` (+ `#45726`, `#45746`) — spherical Bessel asymptotics and the trig reduction — two changes

The report says `j_n(x)` diverges instead of decaying like `1/x`, turning
bad at `x ~ 1e18`.  That is real, but **the headline half was already
fixed upstream**, so the tracker entry looks stale:

* `gsl_sf_bessel_j0_e`, `j1_e` and `j2_e` were switched from
  `gsl_sf_sin_e`/`cos_e` to the system `sin`/`cos` by `cd2dd0519` (2013)
  and `bd5b94b47` (2017).  Measured on this tree, `j0(1e20)` is
  `-6.452513e-21 = sin(1e20)/1e20` — correct.  `jl_e` never used the trig
  functions for large `x` (it goes through the asymptotic Bessel).
* What remained was the **Y family** (`y0`, `y1`, `y2`, and `yl`/
  `yl_array` through their seeds), still on `gsl_sf_sin_e`/`cos_e`:
  `y0(1e20)` returned `+2.86e22` with `GSL_SUCCESS`.

Two changes, deliberately separate, both `[upstream]`:

1. **`bessel_y.c` (commit a08ef2f7f).**  Large-`x` branches use the
   system `sin`/`cos`, matching the *j* functions.  Fixes `#36152`'s
   spherical family and `#45726`.
2. **`trig.c` (commit 8a46ec7cf).**  The root cause: `gsl_sf_sin_e`/
   `cos_e` reduce with a three-term pi/4 and lose the angle above
   `~1e16`.  Replaced above `2^21` by the fdlibm Payne-Hanek reduction.
   Fixes `#45746` and removes the reason the Y functions ever needed a
   workaround.

**Why two commits, not one.**  The `trig.c` change reaches `airy.c`,
`sinint.c` and `legendre_*.c`; the narrow `bessel_y.c` change does not.
The *j* functions were fixed one module at a time upstream for the same
reason, and a reviewer accepting the Bessel fix should not be forced to
accept a trig-kernel rewrite with it.  They can be cherry-picked
independently.

**Verification.**  `y0(1e20)` is now `-7.63970404441728225e-21`, matching
`sqrt(pi/2x)Y_{1/2}(x)` from mpmath; seven Y vectors at 1e18–1e22 fail
with `bessel_y.c` reverted.  The reduction was checked against the
platform libm (itself checked against mpmath) on two million random
doubles over the full exponent range: worst absolute error `2.2e-16`, no
disagreement above `1e-13`.  Eight trig vectors fail with `trig.c`
reverted.  `ctest` 56/56.

**Follow-up (commit 5639c380f).**  The note above said
`gsl_sf_angle_restrict_*` was "not done"; that turned out to be the wrong
call.  The functions had the *same* defect — the three-term `2 pi`
reduction — so `angle_restrict_pos_e(1e12)` was off by `8.6e-7` and every
larger argument was worse, silently.  `gsl_sf_clausen_e` and
`gsl_sf_polar_to_rect` both go through the angle restriction, so the
claim that `clausen` is "unchanged" was only true because it inherited the
error.  Commit `5639c380f` routes both functions through the same exact
reduction above `2^21`; `clausen(1e12)` now returns
`-0.93720759242145504459` (true) instead of `-0.93721134493561908`.  The
`GSL_ELOSS` cutoff is still `0.0625/eps`; only the values below it changed.
18 vectors fail with `trig.c` reverted.

**Not done:** `gsl_sf_sinc_e` is fixed for large `x` as a side effect (it
calls `gsl_sf_sin_e`); no new vector there yet.


### `#40755` (+ `#42042`, `#37209`, `#45265`, `#52927`) - the "Group A" Bessel reports

Five Bessel reports were re-examined on the current build (MSVC x64,
Release; values reproduced through the DLL, references from mpmath at
up to 140 digits and scipy).  They had been grouped as "probably already
fixed by the exact-reduction / Bessel work"; only one was, and not for
that reason.

* **`#37209` - `jl_e` NaN for large `l` - already fixed upstream.**
  Commit `441bc40ff` ("partial fix for bug #37209 + test case") is in
  `savannah/master`, and `test_bessel.c:225` is enabled.  `jl(364,
  36.62)` now returns the subnormal `1.1189e-318` with `GSL_SUCCESS`
  instead of a NaN; `jl(149, 1.0)` returns `GSL_EUNDRFLW` with a finite
  value.  The test passes.  No fork change.

* **`#40755` - sporadic `Jn` NaN - fixed (a43fc0055).**  The branch test
  `GSL_ROOT4_DBL_EPSILON * x > (n*n + 1.0)` overflowed the `int` product
  for `n > 46340`; the wrapped negative value made the test true and
  sent `Jn` to the large-argument asymptotic, which returned `nan` with
  `GSL_SUCCESS`.  `Yn` had the same overflow and returned finite but
  wrong values.  Both now compute the product in double and route to the
  Olver asymptotic, matching scipy.  See FORKNEWS.

* **`#42042` - `Jnu` NaN - fixed (e6e34279a).**  At half-integer order
  the `x >= 2` normalization runs through the unnormalized `J_mu`; for
  `nu = 1/2` that is `J_{-1/2}(x) = sqrt(2/(pi x)) cos(x)`, which is
  zero at `x = 3 pi / 2`.  The backward recurrence produced exactly zero
  and the code divided by it.  The fix uses the endpoint `mu = +1/2`
  when `|cos| < |sin|`.  See FORKNEWS.

* **`#45265` - `J0` error underestimated for `x > 4` - not
  reproducible.**  A 128-point sweep over `[4, 1000]` gave a worst
  `true_err / err = 0.57` (at `x = 29.5`), with the value correctly
  rounded against mpmath.  The report is a single data point on 32 bit
  glibc 1.16; recorded as rejected.

* **`#52927` - `j2` make-check failure - not reproducible.**  The
  reporter's mismatch was `1.3e-20` against a `2.8e-22` error bar on
  glibc 2.12; here the same call is correct to `1.1e-22` inside the same
  bar.  The large-`x` `j2` vectors are already disabled under `#if 0`
  (the `#45730` block) precisely because that error estimate is not
  portable, and `j2` calls the system `sin`/`cos`, so the fork's `trig.c`
  reduction never reaches it.  Recorded as rejected.


### `#43256` + `#68312` — the 6j symbol for large angular momenta — fixed

Two defects in the same function, both from evaluating the Racah sum
directly.

**#43256.** `gsl_sf_coupling_6j_e()` formed `delta()` as the product of
three exact factorials divided by a fourth, and the term loop called
`gsl_sf_fact_e()`, which fails at about `171!`.  The report's symbol is
written with actual angular momenta, `(14,16,16;78,62,76)`, but GSL takes
*twice* those values, so the call that overflows is
`gsl_sf_coupling_6j_e(28,32,32,156,124,152)`.  Its exact value, from the
report's own closed form, is `3.60285887003211998154e-4`; the built
library returned `+Inf` with `GSL_EOVRFLW`.

**#68312.** Below the overflow threshold the alternating sum still
cancels.  At the all-equal case `j = 100` the largest term of the Racah
sum is `3.1e14` times the result, so a double-precision direct sum
retains only one or two correct digits; the reporter measured `~1e3`
relative error.  The 9j symbol and the Racah W coefficient are sums of 6j
values and inherited the error.

**Fix.** `gsl_sf_coupling_6j_e()` now uses the three-term recurrence of
Schulten and Gordon in the first angular momentum, adapted from the
public-domain SLATEC routine `DRC6J`: forward and backward recursions are
run and matched in the interior, the sequence is rescaled as needed, and
the result is normalised by
`(2 L4 + 1) sum (2 L1 + 1) {6j}^2 = 1`; the overall sign is fixed by the
phase convention on the last coefficient.  No factorial is formed and no
cancelling sum is taken.  The now-unused `delta()` and `locMin5()` are
deleted.  The 9j and Racah W are unchanged and inherit the stable 6j.

**Validation.** The full suite passes 56/56.  Independent references came
from sympy's exact `wigner_6j`/`wigner_9j` and from a 60-digit evaluation
of the Racah formula, never from GSL.  The recurrence was first
implemented in Python and reproduced the references to a few ulp over
hundreds of random integer and half-integer cases; the C port was then
compared with sympy through ctypes over four thousand more cases, and the
specific large cases were checked to 40 digits.  The one value the report
gives for the 3j, the closed form `(100 100 100; 0 0 0)`, is produced by
the existing edge recursion and is now pinned by a test.

Negative control: with commit `59fc479e2` reverted, the seven new 6j/9j
vectors fail with `val = inf` and `GSL_EOVRFLW` and the suite reports one
test failure; with the fix in place all pass.

One thing the report gets right by accident and is worth writing down:
the arguments in `(14,16,16;78,62,76)` are actual angular momenta, so the
literal GSL call with those integers is a *different*, perfectly finite
symbol (`-4.3149240273135e-3`); the overflowing call is the doubled one.
The regression vector uses the doubled form and the exact value.

See FORKNEWS; commit `59fc479e2`.


## Group F — non-termination: #31362, #48702, #50459, #21836, #45925

Reviewed 2026-10-03 against the built `build-cmake/gsl.dll` with scratch
programs under `%TEMP%`.  Five reports were grouped as "elliptic/gamma
loops and stalls".  Only one is a genuine non-termination.

### #31362 / #48702 — complete elliptic integrals on NaN

`gsl_sf_ellint_Kcomp_e(GSL_NAN)` returns `GSL_EMAXITER` here in 0.0 ms; it
does **not** loop.  Upstream `9a7cf12e0` already added the `nmax` counters
to the Carlson RC/RD/RF/RJ loops, turning the original infinite loop into a
MAXITER error, and left the `FIXME` in `test_sf.c` noting the code should
be a domain error.  `9b3cce337` is the commit that introduced the NaN test
vectors ("infinite loop for ellint with NAN argument").  #48702 is the same
defect as a 2016 OS-dependent "stall".  Fix: reject NaN with `DOMAIN_ERROR`
before iterating, in Kcomp, Ecomp, Dcomp and Pcomp; `gsl_isnan(k)` not
`isfinite`, so `k = +-inf` keeps its existing `k*k >= 1.0` EDOM path.

Negative control: with the guard removed the vectors return `GSL_EMAXITER`
(and before `9a7cf12e0` did not return at all).  This reverses the earlier
`sin_pi`/`cos_pi` policy of letting NaN propagate — justified in FORKNEWS
because these functions have a bounded domain (`|k| < 1`).

### #50459 — incomplete gamma, a < -2^53

Reproduced: `gsl_sf_gamma_inc_e(-1e25, 0.1)` still running after 8 s.  The
report's diagnosis is exactly right: the negative-a recurrence does
`alpha -= 1.0` until `alpha > a`, and at `|alpha| >= 2^53` the decrement
rounds to `alpha`, so the loop never terminates.  Fix: guard `a < -2^53`
and return 0 with `GSL_EUNDRFLW` (Gamma has underflowed).  Negative
control: the vector hangs without the guard.

### #21836 — exact P + Q = 1

Reproduced `P(0.3,1) + Q(0.3,1) - 1 = -1.44e-15`.  Making P and Q use the
same branch everywhere would force one of them to be a subtraction of the
other in a region where it is the small tail, degrading its accuracy
(e.g. `Q(1,10) = e^-10`).  The fix is therefore scoped to where the
subtraction is lossless: `a >= 0.2`, `x < 20`, `x <= a + 2 sqrt(a)`.
Measured `1 - P_series` against Q's own CF across the window: relative
difference below 1e-13 throughout, and the identity is then bit-exact.
Outside it Q keeps its tuned method.  Negative control: the grid reports
1099 failures with the new branch removed, 0 with it.

### #45925 — "flipped" incomplete gamma functions

Not a bug.  `gsl_sf_gamma_inc` is `Gamma(a,x)`; the complement of `P` is
`Q`.  `Gamma(1,3) = Q(1,3) = 0.049787...` is correct.  The beta half of the
report has the same shape (`I_x(a,b)` vs `I_{1-x}(b,a)` are equal only up
to rounding), and is not a defect either.  Documented in both manual trees;
no code change.

Full CTest suite: 56/56 after the changes.


## Group G — 2F1 / hyperg residual: #21835, #30324, #30510, #41837, #50711, #53876, #53905

Reviewed 2026-10-03 against the built `build-cmake/gsl.dll` (values read
through ctypes) with independent references from mpmath at 40-60 digits.
The reports were retrieved through the Wayback Machine, as the tracker was
not reachable directly.  Four of the reports are fixed by code changes, one
is partial (the `c = a+b` case of #21835), one was already fixed and one is
not reproducible.  #30324, deferred here, was taken later once the filter
was widened.

### #50711, #53905, #21835 (x >= 1 half), #39056 (Wolpert) — one defect, fixed

All four are the same dispatch bug in `gsl_sf_hyperg_2F1_e()`.  The
terminating Gauss series is only reached through the
`fabs(a) < 10 && fabs(b) < 10` branch, and the series itself sits behind the
`|x| < 1` domain check, so:

    #50711  2F1(-1,-13,1,0.651439)   GSL_EUNIMPL   (|b| = 13 over the gate)
    #53905  2F1(-10,2,0.5,0.5)       GSL_EUNIMPL   (|a| = 10 over the gate)
    #39056  2F1(-1,-10,1,0.5)        GSL_EUNIMPL   (Wolpert; disabled vector)
    #21835  2F1(-1,-1,-0.5,1.5)      GSL_EDOM      (x > 1)
    #21835  2F1(0,1,1,11)            GSL_EDOM      (x > 1; fermat a = 0)

When `a` or `b` is zero or a negative integer the Gauss series terminates
and 2F1 is a polynomial, valid for every `x`.  The terminating case is now
dispatched to `hyperg_2F1_series()` before the domain check and outside the
magnitude gate, with the negative-integer-`c` cancellation test preserved.
Fixed in commit `2aa89bae8`; ten vectors added and the two #39056 vectors
re-enabled.  Negative control: 11 failures (10 vectors + summary) with the
dispatch reverted; 0 with it.  The `ldnlwm` vector that shared the `#if 0`
block passes either way and was enabled too.

#39056 was recorded as fixed by the earlier e4c4ac326/882c8361d work, but
that only covered the Monajemi case; its Wolpert vector was still disabled
at `test_hyperg.c:688-696`.  The triage row is corrected to name both cases.

### #21835 (c = a + b near x = 1) — still a limitation, documented

The other half of #21835, `2F1(1,13,14,0.999227196008978)`, is the nearby
singularity `c = a + b` at `x = 1`.  On the current build it returns the
right value, `53.4645144...` (mpmath: `53.464514418791908495`), but with
status `GSL_EMAXITER` (11), so the convergence is exhausted while the
partial sum is accurate to ~9 digits.  `doc/specfunc-hyperg.rst` already
documents this region as `GSL_EMAXITER`.  Adding the A&S 15.3.10 `c = a+b`
formula is a separate, larger change and was not attempted; #21835 is
recorded as **partial**.

### #30324 — extend 2F1 to x < -1 — fixed

The report proposes the MathWorld/A&S transformations to cover `x < -1`,
which GSL rejected at `hyperg_2F1.c` (`x < -1.0`).  Taken once the
eligibility filter was widened: `gsl_sf_hyperg_2F1_e()` now uses the two
Pfaff transformations [DLMF 15.8.1 and 15.8.2] to map `x < -1` to
`z = x/(x-1)` in `(1/2, 1)`, and `gsl_sf_hyperg_2F1_conj_e()` does the same
by summing the transformed (non-conjugate) pair with a new complex Gauss
series.  The description is in `FORKNEWS`; commits `b99ec4e7c` and
`cfe9311cd`.

### #30510 — `hyperg_U(a,b,x)` for x < 0 — fixed

When `a` is not an integer and `b` is a positive integer with `x < 0`, the
A&S 13.1.3 reduction used by `hyperg_U_negx()` is degenerate: its two terms
carry `1/Gamma(1-b)` and `1/Gamma(2-b)`, which are poles at `b = 1` and
`b = 2, 3, ...`.  Measured before the fix:

    U(-0.5, 1, -1)   GSL_EDOM          poch(a,-a) = Gamma(0)/Gamma(a)
    U(-0.5, 2, -1)   GSL_EUNIMPL       b >= 2 integer limit unimplemented
    U(-0.5, 3, -1)   GSL_EUNIMPL

The integer-`b` case `b = n+1` is now evaluated with the DLMF 13.2.9 limit
(the reporter used 13.2.41, the `b`-recurrence route; the direct limit is
simpler), with `ln x` replaced by `ln|x|`.  For `x < 0` the principal value
of `U` is complex - it has a branch cut on the negative real axis - and
`ln|x|` selects its real part, the real continuous solution of Kummer's
equation, which is what a real-valued routine must return.  Fixed in commit
`00859f816`; the new path is limited to `|x| <= 100` and `n <= 170`, outside
which the previous `GSL_EDOM`/`GSL_EUNIMPL` behaviour is kept.

Validation: the formula (with the `1/(n! Gamma(a-n))` prefactor, which is
easy to misread) matches `mpmath.hyperu` at 50 digits for `x > 0`, and for
`x < 0` the `ln|x|` value equals the real part of the principal value,
confirmed by integrating Kummer's ODE analytically from `x > 0` to `x < 0`.
Over a 660-point grid the worst ratio of true error to reported error is
0.61.  Negative control: 13 failures (12 vectors + summary) with the branch
removed.

Not addressed: the non-integer-`b` part of the same region,
`U(-0.5, 2.7, -1)`, still returns `NaN` with `GSL_SUCCESS`; that is a
separate silent-NaN defect.

### #41837 — "bugs in gsl_sf_hyperg_U" — not reproducible on the current build

The three values named in the report are all correct here, against mpmath:

    2F1-style call                     GSL (current)        mpmath
    U(1, -6.67, 1)                     0.1136908226796134   0.1136908226796136
    U(1, -500, 0.1)                    0.0019956088624270   0.0019956088624194
    U(1, -500, 40)                     0.0018481758515476   0.0018481758515484

The follow-up comment's concern - that `hyperg_U_series()` continues to the
infinite sum after a successful finite sum when `1+a-b` is a negative
integer, and then fails inside `gamma.c` - does not reproduce either; a
grid of such calls (`1+a-b = -1, -2`, `beps != 0`, `x > 0`) all return the
correct values with `GSL_SUCCESS`.  The defects were evidently repaired
upstream since the 2014 report, so no fork change is made; the row is marked
fixed (inherited).

### #54998 — confirmed fixed, no further work

The report's `2F1(-0.25, 0.25, 1, 0.9) = 0.920589251382092768` and the
error-estimate case `2F1(-0.25, 0.25, 1, 0.25)` are both covered by the
vectors added with e4c4ac326/882c8361d and pass; the negative-integer
Gamma-prefactor sign defect described in the report's 2025 comments is the
one already fixed there.  Nothing further.

### #53876 — missing `x^(1-c)` in the 2F1 renorm functions — fixed

`gsl_sf_hyperg_2F1_renorm_e()` and `gsl_sf_hyperg_2F1_conj_renorm_e()`
dropped the `x^(1-c)` factor of [A&S 15.1.2] in the negative-integer-`c`
branch, so the result was too large by `x^(c-1)`:

    renorm(1,2,-3,0.4)   2572.016    vs    65.8436
    renorm(5,5,-1,0.5)   4949760     vs    1237440      (4x)
    renorm(5,5,-10,0.5)  1.39e20     vs    6.81e16      (2048x)

Fixed in commit `46b7c412e` with `gsl_sf_pow_int_e()`; the early-termination
branch returns zero and is unaffected.  **The pre-existing test vectors at
`c = -1, -10, -100` were themselves generated from the buggy code and
encoded the missing factor**; they are corrected against the A&S limit
evaluated independently with mpmath (mpmath's `hyp2f1` cannot be used
directly because `Gamma(c)` is a pole, so the ratio of Gamma functions, the
`x^(1-c)` power and the remaining `2F1` are formed separately).  Negative
control: 8 failures (7 vectors + summary) with the factor reverted.

Files: `specfunc/hyperg_2F1.c`, `specfunc/test_hyperg.c`.  Commits
`2aa89bae8`, `46b7c412e`.  Full CTest suite: 56/56.


## Group E — the 2020 specfunc NaN/limit cluster: #55687, #58031, #58060–#58065

All filed by Jackson Vanover, 2019–2020, none with a patch attached, so
every item is a from-scratch reproduction.  The reports are grouped
because they all concern the same question: what should a special
function do with an argument that is non-finite, or a parameter at which
the general formula is singular.

The rule the fork applies, settled by `#57978` and `#57979`:

* a **NaN** argument propagates; it is never a domain error;
* a non-finite argument **with a finite limit** is in the domain
  (`GSL_SUCCESS`);
* a non-finite argument with no limit, or non-physical for that slot,
  is a domain violation (`GSL_EDOM`);
* a removing value (the parameter does not affect the result) is
  returned, which is the same limit semantics the fork already gave
  `1F1_e` in `#58032`.

Reproduction used the built `gsl.dll` under ctypes, so each symptom was
observed rather than inferred.

**Fixed:**

* `#55687` — `gsl_sf_hyperg_1F1_e(1, NaN, -1)` **crashed** with a stack
  overflow.  With `b = NaN` nothing classifies, the `b < 0` branch calls
  `hyperg_1F1_ab_neg()` → `hyperg_1F1_U()`, which re-enters
  `gsl_sf_hyperg_1F1_e()` with the same `b`, and the mutual recursion has
  no base case.  A NaN guard returns NaN/NaN, `GSL_SUCCESS`.  Negative
  control: the test program aborts with `STATUS_STACK_OVERFLOW`
  (`0xC00000FD`) without the guard.  Commit `99a73dd36`.

* `#58064` — only the Si/Ci half is a code defect.  `Si(±inf)` and
  `Ci(+inf)` returned NaN because the asymptotics evaluate `cos/sin(inf)`;
  all three have finite limits (`±pi/2`, `0`) and are now returned.
  `Ci(-inf)` stays `GSL_EDOM` (no limit, and `Ci` is undefined for
  `x < 0`).  The other functions in the report (`gamma_inc_P/Q`,
  `gamma_inc`, `erf`, `erfc`, `poch`) already return the correct finite
  limits; they are documented, not changed.  Commit `3afec212b`.

* `#58065` — `C_n^(lambda)(x)` is identically 0 at `lambda = 0` and
  `n >= 1` (scipy and mpmath agree for `n = 1..5`).  GSL returned the
  limit `2 T_n(x)/n`, including from explicit `lambda == 0` branches in
  `_1_e`/`_2_e`/`_3_e`, and the recurrence in `_n_e` seeded from them.
  All four, and the array form, now return 0.  Three existing vectors
  encoded the old limit values and are corrected.  Commit `01d7e7609`.

**Rejected as not-a-bug:** `#58060`, `#58061`, `#58062`, `#58063` and the
`bessel_zero_Jnu` case, plus `#58031`.  These all report that a NaN is
"dropped" where the function returns a value independent of that
argument (`laguerre_n(0,a,x) = 1`, `U_int(0,1,x) = 1`,
`1F1_int(0,1,x) = 1`, `hermite(0,x) = 1`, `gegenpoly_n(0,lambda,x) = 1`).
Under the fork's rule that is correct propagation: the answer does not
depend on the unknown, so inventing a NaN would be wrong.  `#58031`
(`Kn_scaled(0,+inf) = 0`) is the finite limit `e^x K_0(x) -> 0`, also
correct.  The one real gap was documentation, addressed in `fc4059b12`,
which also adds the NaN rule to the relevant entries and a
`Kn_scaled(0,+inf)` vector.

Recorded in `FORKNEWS` under `[upstream]` (four changes) and `[rejected]`
(one entry), commit `91bdc2d04`.  Full CTest suite 56/56 on MSVC.


## Group I — documentation: #30583, #35032, #36578, #40196, #41605,
## #45234, #58068, #58069, #59912, #66894, #68073, #68592

Reviewed 2026-10-03.  Twelve reports, all in category Documentation.
Two needed no change and are recorded here so they are not
re-investigated:

* `#41605` - the `gsl_histogram_pdf` section is complete in both
  manuals (`histogram.rst` / `histogram.texi`); the online gap the
  report described is gone.
* `#45234` - `specfunc-mathieu.{rst,texi}` already use the SF API
  (`_e` variants, `int` return types).

The other ten were fixed, one commit each, followed by the common
`FORKNEWS` commit `10e4f4df2`:

```
#30583  b0eec8bc6  doc: relate the Legendre and Carlson elliptic forms
#35032  8ef88af6e  doc: gsl_test support functions in usage
#36578  060c2472c  doc: correct the gsl_ieee_env_setup description
#40196  6feb963d3  doc: state the gsl_integration_qag key coercion
#58068  42f33104b  doc: document the gsl_sf_bessel_Jnu domain
#58069  e90d4f002  doc: allow a = 0 in gsl_sf_gamma_inc_Q
#59912  aa796e20b  doc: note the permutation application header files
#66894  dbe316d49  doc: note that randist does not validate arguments
#68073  521bed61d  statistics: correct the gsl_stats_select comment
#68592  20eda0316  doc: document the binary search tree module
```

Method notes:

* Only `#68073` touches a compiled file, and only a comment, so there
  is no behaviour to test with a negative control.  Each documented
  fact was instead reproduced with the built `gsl.dll`: `Jnu_e(1,-1)`
  and `Jnu_e(0,0)` return `GSL_EDOM`; `gamma_inc_Q_e(0,2.2)` returns 0
  with `GSL_SUCCESS` while `gamma_inc_P_e(0,2.2)` returns `GSL_EDOM`;
  `gsl_stats_select` returns the k-th smallest; and the elliptic
  identity was checked with mpmath.
* `#58069` has a subtlety worth keeping: `Q` accepts `a = 0` but `P`
  does not, despite `P = 1 - Q`.  The doc was widened only where the
  code allows it.
* The Sphinx manual builds with the same 5 pre-existing warnings and no
  new ones.  `makeinfo` is not installed here, so the Texinfo edits were
  checked structurally only (balanced directives, consistent node
  chain); they should be rebuilt with GNU Texinfo before being trusted.
* `#68592` and `#35032` carry an explicit note marking the new text as
  an unreviewed, AI-assisted draft, per instruction.
* `bst_test` passes 186476/186476.


## Group H - platform-specific test failures: #39152, #46593, #47028,
## #48915, #49518, #49697, #52127, #52322, #54919, #56843, #59759,
## #67446, #67447, #67705

Reviewed 2026-10-03 as a batch.  None of the fourteen reproduces on any
host the fork builds and tests on, and none carries a code change.  They
share a shape: a make-check or unit-test failure reported against one
toolchain, ABI or architecture, with no demonstrated defect in the
library.

**Hosts.**  MSVC 14.32 x64 (`build-cmake`, Release) and MinGW gcc 15.2.0
x64 (`build-mingw`, Release, static), the latter also with `-O2 -mavx`
(`build-mingw-avx`).  The machine has AVX2 and FMA but no AVX-512.  icc,
ppc64le, AIX and 32 bit builds are not available, so the reports that
need them could not be exercised.  On MSVC x64 every module named below
passes; MinGW passes the same set except `specfunc`, whose failure is a
property of the MinGW C library and not of these reports - see the
method note at the end.

### #39152 - make check errors with Intel icc 13.0.1

The attached logs cover `specfunc`, `ode-initval2`, `poly` and
`multiroots`.  icc is not installed; all four modules pass on MSVC x64
and gcc 15.2 x64.  A 2013 icc failure with no compiler to reproduce on is
indistinguishable from the over-tight test tolerances that have since
been relaxed (e.g. `5432a392c`), so the report is closed as
toolchain-specific.  The `.log` attachments were not downloaded by the
inventory script, so only the module list is available.

### #46593, #47028, #52322 - multifit on 32 bit, ppc64le and 32 vs 64 bit

`multifit.test` passes on MSVC x64 and gcc 15.2 x64.  No 32 bit or
ppc64le host is available.  #52127 below is the same class and was
checked against the reporter's own program; the differences these reports
describe are at the last-ulp level between ABIs, not a wrong value.
#52322's attachments were not downloaded and the report names no incorrect
number, only a difference.

### #48915 - some test failures on AIX (GSL 2.1.91)

Names `splinalg`, `linalg`, `multilarge_nlinear`, `rng`, `poly` and
`specfunc`, all with logs and no patch.  No AIX host.  The five
non-`specfunc` modules pass on MSVC x64 and gcc 15.2 x64; `specfunc`
passes on MSVC and fails on MinGW only for the libm reason below.  The
report is from 2016 and predates the tolerance work; closed as
platform-specific with no reproduction.

### #49518 - "bug in matrix/vector tests" under MSVC x64

**Already fixed; no fork change.**  The report is real and specific: five
test files used

    char filename[] = "test.XXXXXX";
    #if !defined(_WIN32) ... mkstemp ... #else char *fd = _mktemp(filename); ...

and on MSVC x64 `_mktemp()` was called with no prototype (its declaration
lives in `<io.h>`), so the returned pointer was truncated through an
implicit `int` return and the test crashed.  The current tree contains no
`_mktemp()` or `mkstemp()` anywhere: `matrix/test_source.c`,
`vector/test_source.c`, their complex variants and `spmatrix/test_source.c`
now use a fixed `test.dat` / `test_static.dat` filename and `unlink()`.
`matrix`, `matrix.test_static`, `vector`, `vector.test_static` and
`spmatrix` all pass on MSVC x64.  There is nothing left to fix.

### #49697 - linalg test fails with gcc and -mavx

`linalg.test` passes on MinGW gcc 15.2 x64 built with `-O2 -mavx`.  The
host exposes AVX2 and FMA but no AVX-512, which is the widest case the
report could need.  The 2016 failure is almost certainly the unscaled
Hilbert `gsl_linalg_cholesky_invert` check, whose tolerance the fork has
already raised from `256` to `512` times `N * GSL_DBL_EPSILON` for
Savannah #67445 (`linalg/test_cholesky.c:241`).  Recorded as not
reproducible / already covered.

### #52127 - gsl_eigen_nonsymm differs between 32 and 64 bit

The attached `nonsymm.c` was built against the fork's library and run.  On
x64 it prints

    Real                 Imaginary
    3FF0000000000000     40158A68A4A8D9F4
    3FF0000000000000     C0158A68A4A8D9F4
    3FF0000000000000     40158A68A4A8D9F5
    3FF0000000000000     C0158A68A4A8D9F5

i.e. `1.0000000000000000 +- i 5.3851648071345`, which is exactly the
result the report predicts for 64 bit.  The report's 32 bit values
(`1.00000000000000090`, `0.99999999999999956`) are the outlier, off by
one or two ulp of the real part.  The matrix is `1 + 2i + 3j + 4k`; its
eigenvalues are `1 +- i sqrt(29)`, and the x64 answer agrees with that to
an ulp.  Bit-identical eigenvalues across ABIs are not a GSL guarantee.
The second attachment is an unrelated `graph_sampler` configuration file
and adds nothing.

### #54919 - gsl 2.5+ test fails with icc (2016.4 and later)

No icc.  The named modules pass on MSVC x64 and gcc 15.2 x64.  Same
disposition as #39152.

### #56843 - linalg eigen fail on non-x86 hardware (accuracy)

`eigen.test` passes on MSVC x64 and gcc 15.2 x64.  The report is about
accuracy differences on "non-x86" hardware; the only such target the fork
builds is arm64 macOS, which is covered by CI.  Differences at the
`1e-16` level are expected when the intermediate precision differs (x87
on x86-64 against rounded `double` on arm64) and are not a defect.

### #59759 - spmatrix test fails on x86_64

The report body is about 32 bit machines.  `spmatrix.test` passes on
MSVC x64 and gcc 15.2 x64.  No 32 bit host; closed as ABI-specific.

### #67446, #67447 - multilarge_nlinear and spmatrix under gcc 14.2.1

Both pass under MinGW gcc 15.2 x64 (and MSVC x64).  These are the most
recent of the batch (August 2025) and were the most likely to reproduce,
but the exact gcc 14.2.1 build could not be reconstructed and gcc 15.2
does not show either failure.

### #67705 - ttest failure in linalg/QR_solve_r random

`linalg.test` passes on MSVC x64 and gcc 15.2 x64.  The `QR_solve_r
random` vectors draw from `gsl_rng_default`, whose default seed is fixed,
so the input matrices are deterministic; the run is not flaky here.  No
failing `N` could be produced, so the report is left as not reproducible.

### Method note - MinGW's libm, not a Group H report

On MinGW gcc the `specfunc` suite fails in the spherical Bessel Y
functions: the fork-added vectors at `1e18`-`1e22` (Savannah #36152) are
wrong by several percent, and the upstream vectors at `2^32` trip the
value/error consistency check.  The cause is the MinGW C library, whose
`sin`/`cos` argument reduction is already wrong at `1e18`:

    x        MinGW sin            accurate sin
    2^32     -0.46198657951383493 -0.4619865795138349
    1e18     -0.99281610405300347 -0.9929693207404051
    1e20     -0.74692189125949293 -0.6452512852657808
    1e22      0.46261304076460175 -0.8522008497671888

MSVC's UCRT and glibc reduce these correctly.  The Y functions delegate
to the system `sin`/`cos` for large `x` (by design, matching the *j*
functions), so they inherit the MinGW defect.  MinGW is not one of the
fork's build targets (MSVC, Linux, macOS), so this is recorded but not
acted on here.


## Group J - distribution inverse / quantile bugs: #42502, #45924, #67058

Reviewed 2026-10-03.  Three reports of the same shape as the applied
`#39057` (gamma inverse): a quantile or aggregate that returns a wrong
value or `NaN` for a corner of its parameter space.  Two were real
defects, fixed; one was a false report.  All three were reproduced
against the built `gsl.dll` through ctypes before anything was changed,
and only reads of the built library were used as references - scipy
1.15 for the beta inverses and the R/numpy convention for the empty
aggregate.

### `#45924` - `gsl_cdf_beta_Pinv` "inverse failed to converge" - fixed

All ten cases in the report (`a = 8000, b = 2000`; `a = 630, b = 9370`;
`a = 5000, b = 5000` at `P = 0.005, 0.995, 5e-5, 0.99995`) returned
`NaN` on the built DLL.  scipy agrees with the Boost values quoted in
the report.

Three defects, one root cause.  `gsl_cdf_beta_Pinv` bisected to an
absolute tolerance in `P` (`xtol = Ptol = 0.01`) and then ran an
unsafeguarded Newton iteration.  For large `a,b` the initial small-x
approximation already satisfies `|Px - P| < 0.01` while being nowhere
near the root (for `a = 8000, b = 2000`, `|Px-P| = 0.005` on the first
probe), so `bisect` returned it immediately and the Newton iteration
could not recover.  Second, for `a << 1` the quantile lies hundreds of
decades below the mean; a bisection on `x` cannot cross that in 64
steps.  Third, `gsl_cdf_beta_Qinv` formed the upper tail as
`1 - gsl_cdf_beta_Pinv(Q, b, a)`, which loses every digit by
cancellation when the result is small - measured relative errors up to
`3e292` on the random grid.

The solver is reworked around `t = logit(x)`, the analogue of the
`log x` used for the gamma inverse.  In both the lower-tail
(`F = P(x) - target`) and upper-tail (`F = target - Q(x)`) forms `F` is
monotone increasing in `t`, so `[log(DBL_MIN), -log(DBL_MIN)]` is
always a bracket.  A Newton step from the density is accepted when it
stays inside the bracket and changes `x`; otherwise the bracket is
halved.  Convergence is accepted when the bracket closes, when the
Newton step is below the resolution of `t`, or when `x` can no longer
change (many `t` map to the same double `x` near 0 and 1).  A root
below `GSL_DBL_MIN` returns zero, matching `gammainv.c`.  Both entry
points now invert whichever tail is no larger than one half, so `Qinv`
never subtracts from one.

The iteration needed four tries to get right, all caught by the random
grid rather than by reading: the first `x`-based tolerance stalled when
a stale bracket endpoint kept `xhi - xlo` above the threshold; a pure
`t`-width tolerance stalled when `x` was pinned but `t` was not; an
`xnew == x` break fired prematurely in the middle of the search; and a
small Newton step with a stale endpoint needed an explicit acceptance
test.  The final rule accepts any of the three conditions.

Verification: 73085 `Pinv` and 99658 `Qinv` random points (`a,b` in
`[1e-4, 1e4]`, tail probability from `1e-12` to `1/2`) have no
non-finite results and a worst relative error of `7.4e-11`, against
`TEST_TOL6 = 2.3e-10`.  The ten report cases agree with scipy to about
`4e-14`.  `gsl_cdf_fdist_Pinv/Qinv` delegate here and inherit the fix;
two large-`nu` fdist vectors are added.  Negative control: with
`betainv.c` reverted, `cdf_test` aborts at the first new vector
("inverse failed to converge", exit `0xc0000409`), and with the error
handler off all ten report cases return `NaN`.  Full suite 56/56.
Files: `cdf/betainv.c`, `cdf/test.c`.  Commit `51a63cbd5`.

### `#67058` - empty `gsl_stats_mean` / `gsl_stats_sd` return 0.0 - fixed

Measured: `gsl_stats_mean`, `gsl_stats_variance` and `gsl_stats_sd` of
a zero-length array all returned `0.0`.  The manual defines the mean as
`(1/N) sum x_i`, undefined at `N = 0`; R and numpy return `NaN`.  The
cause is the recurrence `mean += (x_i - mean)/(i+1)`, which leaves the
accumulator at zero when the loop never runs; the variance family
inherits the zero and the `1/(N-1)` factor cannot correct it.

`gsl_stats_mean` and the shared `compute_variance` helper (so
`variance_m`, `sd_m`, the fixed-mean pair, and `variance`/`sd` through
them) now raise `GSL_EBADLEN` for `size == 0`, returning `NaN` with the
error handler disabled.  `gsl_stats_tss` is deliberately unchanged: the
empty sum is genuinely zero.  The behaviour is recorded in both
manuals.  Negative control: with the two guards reverted,
`statistics_test` reports 99 new failures ("... empty"), all in the
empty-array checks; with them applied the suite is clean.  Files:
`statistics/mean_source.c`, `statistics/variance_source.c`,
`statistics/mean.c`, `statistics/variance.c`, `statistics/test.c`,
`statistics/test_float_source.c`, `statistics/test_int_source.c`,
`doc/statistics.rst`, `doc_texinfo/statistics.texi`.  Commit
`68fc3c752`.

### `#42502` - `gsl_cdf_ugaussian_Pinv(0.5)` "returns 1" - not a bug

The report's program includes only `<gsl/gsl_math.h>` and `<stdio.h>`,
not `<gsl/gsl_cdf.h>`, so `gsl_cdf_ugaussian_Pinv` is implicitly
declared as returning `int`.  The function does return the double `0.0`
(it is in `xmm0`); the `1.000000` is the second `%f` reading a stale
vararg slot.  On the built DLL `Pinv(0.5)` is exactly `0.0`, and the
call is already covered by `cdf/test.c:498`
(`TEST (gsl_cdf_ugaussian_Pinv, (0.5), 0.0, TEST_TOL0)`), so no code
change is made.  Recorded in `FORKNEWS` as `[rejected]`.

Recorded in `FORKNEWS` under `[upstream]` (two changes) and `[rejected]`
(one entry).  Full CTest suite 56/56 on MSVC x64.


## Group K - special-function accuracy, 2025-2026: #68479, #68283, #43259, #66993, #67494

Reviewed 2026-10-03 on branch `group-k` (base `b999b7e05`), in a
separate worktree from the parallel Group L work.  Every reproduction
is against the built DLL; references are mpmath at 50-60 digits.  The
applied changes are recorded in `FORKNEWS`.

### #68479 - gsl_ran_binomial, BINV seed - fixed

The small-mean branch seeds `f(0) = (1-p)^n` with `gsl_pow_uint(1-p,n)`.
The rounding of `q = 1-p` is amplified by `n`, so the seed's relative
error is `~n*eps`:

    p=1e-9  n=1e9  q^n = 0.36787945073036543   exp(n log1p(-p)) = 0.36787944098750258
    p=1e-8  n=1e9  q^n = 4.5399925828314055e-05 exp(n log1p(-p)) = 4.5399927492488394e-05

`exp(n*log1p(-p))` agrees with mpmath; the seed is the base of every
probability in the inverse-CDF walk, so the whole distribution was
perturbed.  Fixed as posted (commit 7376dc316), with a regression test
that returns a fixed uniform in the gap between the two seeds.

### #43259 - accuracy problems in specfunc - partly fixed

Fifteen inputs across nine functions.  The report's "actual error"
usually equals the magnitude of the GSL value, i.e. the reporter's
reference underflowed to 0, so the figures are not evidence; each case
was re-derived here.

| function | current build | reference | verdict |
|---|---|---|---|
| `clausen` x2 | `0.7414877353776` / `0.7813604462495` | match to ~4e-16 | already fixed (5639c380f) |
| `zeta`, `eta` | match to ~4e-16 | err now >= actual | already fixed |
| `psi_1` x2 | `GSL_EDOM` | arguments are poles (negative even integers) | correct |
| `exprel_2` | rel err `8.9e-14`, reported `4.4e-16` | mpmath | **fixed** (c2ba2e271) |
| `gamma_inc_Q` | rel err `1.6e-5` | `6.3191220136968329e-11` | **fixed** (18652f8c1) |
| `hyperg_0F1` x3 | `1.15e198`, `1.14e168`, `2.23e157` | `1`, `2`, `1.7893298094111404e34` | **fixed** (5fafddc63) |
| `pochrel` x3 | `~203`, `~208`, `~133` | `a ~ -1e89` sits on a `Gamma` pole | degenerate |
| `ellint_P` | `4.148e-119` | `n ~ 1.3e251` puts a pole in `(0,phi]` | degenerate |

The `hyperg_0F1` failure is the instructive one: for `c < 1` the code
used `sin(pi(1-c))`, and with `c ~ 1e-215`, `1-c` rounds to `1`, so
`sin(pi)` is `1.2e-16` instead of `sin(pi c) ~ 1e-214`; the `K` term is
then wrong by ~120 orders of magnitude.  The fix passes `c` to the
helpers and uses `gsl_sf_sin_pi`/`gsl_sf_cos_pi`.  `exprel_2` and
`gamma_inc_Q` are cancellation, not mis-classification: the direct
`exp(x)-1-x` and the `1-P` complement lose digits for small arguments,
and both now use the cancellation-free series when appropriate.

A 475-point `gamma_inc_Q` grid (`a` 1e-1..1e-19, `x` 1e0..1e-24) has no
point worse than before with the new branch threshold `a|log x| < 0.1`,
and `P + Q = 1` is still exact in its window.

### #68283 - gsl_rstat_skew/kurtosis "correction" - rejected

The proposed guards (`n>2 && M2>0`, `n>3 && M2>0`) are inconsistent
with the module's own contract: `gsl_rstat_*` is the online equivalent
of `gsl_stats_skew`/`gsl_stats_kurtosis`, and `rstat_test` asserts that
equivalence.  Measured on the built DLL:

    n=1 [5]         skew=nan  kurtosis=nan
    n=2 [1,2]       skew=0    kurtosis=-2.75
    n=3 [1,2,3]     skew=0    kurtosis=-2.333
    n=4 constant    skew=nan  kurtosis=nan

`gsl_stats_skew` returns `NaN` for `n=1` and `gsl_stats_kurtosis` the
listed values for `n=2,3`.  The proposed change would fail test1 at
`j=1` (0 vs `NaN`) and `j=2,3` (kurtosis 0 vs -2.75/-2.333), because
`gsl_test_rel` treats value-vs-`NaN` as a failure.  No code change.

### #66993 - pow_int issues - rejected (one doc note)

Four claims, none reproducing on the current build:

    (0.5)^-9        = 512.0 exactly   (report: 512.0000000000010232)
    (-0.3)^1        = -0.3 exactly    (report: -0.3000000000000002)
    (-0.0)^-1       = -inf            (report: -7205759403792794)
    0.0^-1          = +inf, GSL_EOVRFLW (report wanted GSL_EZERODIV)

The deprecation/rationalisation proposal is an API change and out of
scope.  The only valid point is that `doc/specfunc-pow-int.*` says the
functions never check overflow, while `specfunc/pow_int.c:42` returns
`GSL_EOVRFLW` for a zero base with a negative exponent; that exception
is now documented (afd7f7d6c).  Changing the code to `GSL_EZERODIV`
would be a behaviour change with no demonstrated defect.

### #67494 - consistency of pow_int usage - rejected

`gsl_pow_int(2.0, m)` and `gsl_sf_pow_int(2.0, m)` in `hermite.c:232`
return identical values and no error estimate is used; purely
stylistic.  No change.

Full CTest suite 56/56 on MSVC x64.  Negative controls for each applied
change are in `FORKNEWS`.


## Group L - non-termination / domain handling: #21837, #31426, #40176, #52351, #59911, #59913

Reviewed 2026-10-03.  Six reports with a non-termination, array-index,
domain or error-code theme.  Four were real and are fixed or corrected;
one was already fixed upstream; one is a property of the caller's
integrand.  All were reproduced against the built DLL before any
change, in a `group-l` worktree branched from `b999b7e05`.

### `#52351` - "akima.c array indexing" - real, fixed

There is no out-of-bounds access.  The defect is the handling of the
Akima weight denominator in `akima_calc()`: when the two neighbouring
differences are equal the weight is indeterminate and GSL substituted
one of the adjacent slopes (`m[i]` for both tangents), which makes the
spline non-C1.  The correct limit is the average of the two adjacent
slopes.  The reporter's proposed `m[i+1]` is also wrong: it removes the
corner but overshoots the minimum.

Method: transcribed the routine into Python and compared with the
independent scipy `Akima1DInterpolator(method="akima")`.  The fixed
code matches to `4.4e-16` on the demo2 data; the old code and the
reporter's variant are off by `0.140625`.  Negative control: with
`akima.c` reverted the six new vectors in `interpolation/test.c` all
fail (`-1.5` vs `-1.625` at `x = 4.5`).  Full suite 56/56.
Files: `interpolation/akima.c`, `interpolation/test.c`.  Commit
`0efb5f87f`.

### `#31426` - infinite loop in `gsl_eigen_symm` - real, fixed

The 1024x1024 `mat.bin` attachment is no longer reachable, but the
failure class is reproducible: a tridiagonal matrix whose entries are
subnormal makes the chopping threshold
`GSL_DBL_EPSILON*(|d_i|+|d_{i+1}|)` underflow to zero, so
`chop_small_elements()` never removes the subdiagonal, the QR sweep
reaches an exact fixed point and `b` never decrements.  For `N = 16`
with all entries `GSL_DBL_MIN`, `d` and `sd` are unchanged after the
first sweep.

Fixed by rescaling the referenced triangle by a power of two (exact;
the Householder and QR steps are scale-equivariant) and adding a
per-block sweep bound that returns `GSL_EMAXITER`.  Verified with a
faithful Python transcription, `numpy.linalg.eigvalsh` and the closed
form `1 + 2 cos(k pi/(N+1))`.  Negative control: with
`symm.c`/`symmv.c` reverted `eigen_test` hangs (12 s timeout); with
them it passes, and the standalone reproducer returns `GSL_SUCCESS`.
`GSL_EMAXITER` is documented in both manuals.  Full suite 56/56.
Files: `eigen/symm.c`, `eigen/symmv.c`, `eigen/test.c`,
`doc/eigen.rst`, `doc_texinfo/eigen.texi`.  Commit `7d586d89b`.

### `#21837` - `solve_symm_tridiag` on a zero diagonal - documented, request rejected

Reproduced the report with a returning handler: the zero-diagonal
example returns `GSL_EZERODIV` (12), not a solution.  The maintainer's
2007 reply stands: a solution needs a permutation (`slatec/dgtsl.f`),
a new algorithm and out of scope.  The report did expose a
documentation defect - both manuals claimed `GSL_ESING` while the code
returns `GSL_EZERODIV` - which is corrected.  Files: `doc/linalg.rst`,
`doc_texinfo/linalg.texi`.  Commit `9418afa5f`; the permutation request
is recorded `[rejected]`.

### `#40176` - "possible error in poly test suite" - already fixed

The 15th-order polynomial failure (conjugate pair interchanged) is
#39055.  `poly/test.c` is identical to `savannah/master`, which
contains `9cc12d037` (sort by real then imaginary part, negative
imaginary root listed first) and `0466df866` (relative check).  The
vector passes on the build; nothing to do.  Recorded `[rejected]`.

### `#59911` - qagiu on `cosh(x) exp(-cosh(x))` - not a bug

The `qagiu` transform probes `x = (1-t)/t` as `t -> 0`.  For
`x > ~710.5` the caller's `cosh(x)` overflows to `+Inf`, `exp(-...)`
underflows to 0 and the product is `NaN`.  The built library returns
`GSL_EMAXITER` with a NaN result (the report saw "bad integrand
behavior").  Not something the library can rescue; the thread's own
advice (substitute `u = sinh x`) is right.  Recorded `[rejected]`.

### `#59913` - cquad SIGFPE - latent defect fixed

The original program is unavailable and a 2025 follow-up could not
reproduce it, but `cquad.c` computes `ncdiff / nc` unguarded at the
initial error estimate (the later, identical test is guarded by
`nc > 0`).  For `f(x) = 0`, `nc = 0` and the quotient is `0/0`
(`FE_INVALID`), which kills the process if that trap is enabled (MSVC
exit `0xc0000090`).  Added the `nc > 0` guard and a trap-based
regression test.  Negative control: with `cquad.c` reverted
`integration_test` is killed by the trap; with it the test passes and
the standalone reproducer returns `GSL_SUCCESS`.  Full suite 56/56.
Files: `integration/cquad.c`, `integration/test.c`.  Commit
`f20496de6`.

Recorded in `FORKNEWS` under `[upstream]` (four changes) and
`[rejected]` (three entries).  Full CTest suite 56/56 on MSVC x64.


## Group M — interpolation / histogram / misc correctness: #38548, #68379, #47193, #67301, #68415, #65932, #63519

Seven reports, each reproduced (or refuted) against the built library
before a verdict was recorded.

### `#38548` — histogram rounding for integer data — fixed

The report is real and reproduces exactly: `gsl_histogram_calloc_uniform
(60, 0, 60)` gave `range[31] = 31.000000000000004`, so `find(31.0)`
returned bin 30.  `make_uniform()` built each limit as
`((n-i)/n)*xmin + (i/n)*xmax`, rounding the coefficients first.  The
same helper is duplicated in `init2d.c`; the two-dimensional histogram
had the identical defect (`xrange[31] = 31.000000000000004`).

The posted/mailing-list formula was the midpoint
`0.5*(xmin+i*dx) + 0.5*(xmax-(n-i)*dx)` with `dx=(xmax-xmin)/n`.  It is
exact for the reported case, but comparing three candidates against
exact rational arithmetic showed that it moves the correctly rounded
limits of a 10-bin histogram over `[0,1]` (`0.1`, `0.2`, `0.4`, `0.6`)
by an ulp, as does the single `xmin + i*dx`.  That regression is caught
by the deterministic `test1d_resample` test.  Only
`xmin + (xmax-xmin)*i/n` is exact for both cases, and it is what the
fork uses; `range[0]`/`range[n]` are pinned because the division can
still move them.

Verification: exact-rational reference over integer spans (exact for
every divisible span, monotone in ~27k randomised cases); new vectors
in `test1d.c`/`test2d.c`; negative control with the old formula reports
5 failures (3 one-dimensional, 2 two-dimensional).  `ctest` 56/56.
Commit `f9c977f52`.

### `#68379` — histogram accessor index type — fixed

`gsl_histogram_max()` and the 2D `xmax()`/`ymax()` narrowed `h->n`
(`size_t`) to `int` before indexing.  Latent for any allocatable
histogram, but corrected to `size_t`.  Commit `ae60e2883`.

### `#47193` — `gsl_ran_poisson_pdf` with `mu=0` — already fixed

Duplicate of #43326, fixed earlier in this fork (`4f9f4f4fc`).
Verified: `gsl_ran_poisson_pdf(0,0)=1`, `(1,0)=0`.  The earlier fix had
no test, so one was added; with the `mu==0` guard disabled it is the
only failure in `randist_test`.  Commit `5346bba43`.

### `#67301` — the second `find()` on `r2` — rejected

The report asks for `find(p->nx*p->ny, p->sum, r2, &k)` after the
existing lookup on `r1`.  That would be a bug: `find()` locates a value
in the cumulative distribution, whereas `r2` is the within-bin `y`
fraction.  A second lookup would overwrite the bin selected by `r1` and
combine two different bins.  Measured on the build, `r2 = 1.5` simply
returns a `y` outside the selected bin; it does not crash.  Only the
function description was clarified (commit `2b15efbff`).

### `#68415` — variance at `n = 1` — rejected

The report's premises do not reproduce.  It claims
`gsl_stats_variance_m()` "returns 0.0" and `gsl_stats_pvariance()`
"returns `GSL_EDOM`", but on this build all of variance, covariance and
pooled variance return `NaN` for `n = 1`: the sum of squares is exactly
zero and is multiplied by `n/(n-1) = 1/0`.  The three functions are
already consistent, and `NaN` is the correct result for an estimator
that is undefined with one observation.  Both manuals now state the
`n >= 2` requirement (commit `65bc71664`).

### `#65932` — complex tridiagonal solvers — fixed

A duplicate of #60457, a feature request for complex tridiagonal
solvers.  Both were taken once the filter admitted new public API; the
four solvers are in `linalg/tridiagcomplex.c` (commit `dd5e0a0f2`), with
numpy reference solutions and residual checks.  See `FORKNEWS`.

### `#63519` — `gsl_root_fsolver_set` straddle error — fixed

The request is to return `GSL_EINVAL` without invoking the error
handler.  `GSL_ERROR` already returns the code - under
`gsl_set_error_handler_off()` the call returns `GSL_EINVAL` (4),
confirmed against the DLL - but it also calls the handler, whose
default action is to abort.  That uniform contract is by design and is
documented in both manuals (commit `aaecb86e6`).  The per-call variant
proposed as `gsl_root_fsolver_set_with_values()` (#66576) was taken
later, once the eligibility rule was widened; it avoids re-evaluating
an expensive `f` and is implemented without an ABI change (`33b79b362`).

Verdicts recorded in `SAVANNAH_TRIAGE.md`; the applied changes and the
rejections are in `FORKNEWS`.  Full CTest suite 56/56 on MSVC x64.


## Group N — build / portability / test-quality: #36197, #39120,
## #39165, #41457, #50382, #55965, #63927, #68518

Eight reports, mostly about the build rather than the numerics.  Each
was checked against the tree before a verdict was recorded.

### `#36197` — reserved identifier violation — fixed

The report is real: 278 tracked headers, `.in` files and two design
manuals used include guards of the form `__GSL_MATH_H__`.  An identifier
that begins with an underscore and an uppercase letter is reserved to the
implementation, and in C++ any identifier containing a double underscore
is reserved too, so a consumer that included a GSL header was relying on
names it does not own.

The tracker discussion (Rhys Ulerich, comment #3) settled on dropping
only the leading `__GSL`; that is what was done, so the guards are now
`GSL_MATH_H__` and keep their trailing underscores.  The private guards
`__ERROR_CBLAS*`, `__ROOTS_H__`, `build.h` and the two `ringbuf.c`
includers were renamed under the same rule, as were `gsl_version.h.in`
and the `contrib/wigner.h` guard.  `__BEGIN_DECLS`/`__END_DECLS` were
deliberately left alone (`gsl_mode_t`/`gsl_prec_t` are unwinnable public
API and were never in scope).

The two attached patches (`36197a.diff` 118 KB, `36197b.diff` 179 KB)
were **not** used: the sweep recorded both as `dirty`, they are thirteen
years old, and the tree has moved a long way since.  The rename was
redone mechanically, replacing only the enumerated guard tokens so that
`__cplusplus`, `__GNUC__` and `__attribute__` are untouched.  The root
`gsl/` directory is a set of build-time symlinks and needs no edit.

Verification: clean CMake rebuild and 56/56 ctest on MSVC x64; every
installed public header compiled twice on its own (251 OK; the 15
`spmatrix` headers are not standalone-includable, pre-existing);
negative control with the old `#define` restored in
`gsl_vector_double.h` gives "conflicting types for 'gsl_vector'".
Commit `e40dc7ceb`.

### `#39120` — possible removal of some files — fixed

The five files that no build references were removed: `err/env.c`,
`err/warn.c`, `rng/g05faf.c`, `poly/norm.c` and `poly/inline.c`.
`err/ChangeLog` already records `warn.c` as removed in 2003, and
`poly/inline.c` duplicates the built `poly/eval.c`.

`ode-initval2/modnewton1.c` must stay: `rk1imp.c`, `rk2imp.c` and
`rk4imp.c` `#include` it.  `matrix/matrix.c` is built.  The report's
rename of `eval.c`/`matrix.c` to `inline.c` was rejected as cosmetic and
hostile to the checked-in source lists.  `statistics/wcovar.c`, also
named, was already absent.  Commit `14920fa04`.

### `#39165` — always-true conditionals — rejected (stale)

Both sites the report names are already gone.  The `gegenbauer.c`
expression was superseded by the #58065 rewrite, and `driver.c` lost its
`hmin >= 0.0 || hmin < 0.0` test in bzr revision 4829.  The `hmax`
guard that remains (`hmax > 0.0 || hmax < 0.0`) is not the same
construct: it is false for `0` and for `NaN`, so it is not always true.

### `#41457` — valgrind errors in `matrix/test.c` — rejected (already fixed)

Fixed upstream on the same day the report was filed, by commit
`6ed874986` (the #39101/#39102 fix), which added
`memset(m->data, 0, MULTIPLICITY * n1 * n2 * sizeof(ATOMIC))` to
`matrix/init_source.c`.  The memset covers the `long double` padding
bytes that `long_double_fwrite` reads.  `valgrind` is not available on
the Windows review machine, so the fixed source was inspected rather
than re-run.

### `#55965` — PCG random number generator — fixed

A new generator is a new algorithm and a new public API, which the
eligibility rule originally excluded.  It was taken once the rule was
widened: `rng/pcg.c` provides PCG32 as `gsl_rng_pcg32` (`3fd2c7475`),
implemented from the published algorithm rather than from the Apache-2.0
reference code.

### `#50382` — CMake and NuGet for Windows — rejected (partly done)

The fork already builds with CMake on Windows, Linux and macOS, which is
the half that matters; NuGet packaging is new distribution surface and is
out of scope.  No code change.

### `#63927` — `gsl-without-cblas.pc` — fixed

`gsl-config --libs-without-cblas` already existed, but pkg-config could
not express it.  A `gsl-without-cblas.pc.in` was added and installed by
both builds (`442b7cbd2` autotools, `9ad4ed051` CMake), and
`pkgconfig.test` now checks that the file does not pull in CBLAS.  The
`Libs:` line is `-L${libdir} -lgsl`.  `pkg-config` is not installed here,
so the scripted check has to run in Linux CI.

### `#68518` — stale `configure.ac` construct — fixed

`DISCARD_POINTER` was a `do { } while` macro that only silenced
unused-parameter warnings.  The three uses (`monte/vegas.c`,
`interpolation/akima.c` twice) are now `(void)` casts and the macro is
gone from `configure.ac` (`46fa69ae1`); the CMake template
`cmake/config.h.cmake` lost it in the separate fork commit
(`2afdbb721`).  Two commented-out references in
`interpolation/steffen.c` are left as comments.  Verified with a CMake
build, 56/56 ctest, and `gcc -Wunused-parameter` (negative control
without the cast reports "unused parameter 'xu'").

Verdicts recorded in `SAVANNAH_TRIAGE.md`; the applied changes and the
rejections are in `FORKNEWS`.  Full CTest suite 56/56 on MSVC x64.


## Group O — root finding and ODE: #39713, #42219, #42220, #50712, #66849, #30540, #30947

Reviewed 2026-10-04 against the built `build-cmake/gsl.dll` (MSVC x64),
with the reporters' own programs compiled against it.  Two reports are
real latent division-by-zero defects and are fixed; one final-time
rounding defect is fixed in both ODE interfaces; one is already fixed
upstream; the three remaining were feature-shaped or too deep to fix and
were recorded as rejected.  Of those, `#30540` was taken later once the
filter was widened.  The fixes are commits `8640b846e`, `dda917296` and
`fa1622e11`; the full `ctest` suite is 56/56.

### `#42219` + `#42220` — division by zero at a root — fixed

Reproduced with the reporters' own programs, compiled against the
built DLL:

    bug_gnewton.c   root = -nan(ind), "the iteration has not converged";
                    GSL_IEEE_MODE=trap-common -> 0xC0000090 (SIGFPE)
    bug_hybrid.c    all four hybrid solvers returned the correct root
                    in the default mode, but raised 0xC0000090 under
                    trap-common

`#42219`: `gnewton_set()` stores `phi = enorm(f)` from the fdf
callback, while `gnewton_iterate()` obtains the next residual from the
plain f callback.  When the initial point is a root of fdf but not of
f, `phi0 = 0`, so the relative step reduction computes
`theta = phi1/phi0 = inf` and the iterate becomes NaN.  The reduction
is now guarded with `phi0 > 0.0`.

`#42220`: at a root the dogleg step is zero, `pnorm = 0`, and the
rank-1 update divides `(qtdf - rdx)/pnorm` and `diag^2*dx/pnorm` by
zero, poisoning the QR factors.  The hybrid iterations now return
`GSL_SUCCESS` immediately when `fnorm == 0`, and `compute_wv()` leaves
the update vectors at zero when `pnorm == 0` (defensive; the early
return already covers it because `Q` is orthogonal, so `qtf = 0` iff
`fnorm = 0`).  `broyden.c` has an adjacent `lambda == 0`
`GSL_EZERODIV`, but that is not part of this report and the reporter
left the broyden case commented out; it is left unchanged.

Regression tests in `multiroots/test.c` run all four hybrid solvers at
the root and gnewton with the mismatched callbacks, with the
invalid-operation trap enabled (the hybrid failure is otherwise hidden
by the residual test).  Negative control: with the four files reverted
the test process is killed by the trap (`0xC0000090`); with the fix the
suite is 77/77.

### `#66849` — `gsl_odeiv2_evolve_apply` past the final time — fixed

Reproduced with a right-hand side that records the time range it is
called over.  Integrating backwards from 1 to `1e-12` with a suggested
step larger than the interval:

    t0 + (t1 - t0) = 9.999778782798785e-13 < t1 = 1e-12
    RHS seen down to 9.999999960041972e-13

The final step is clamped to `h = t1 - t0`, but the stepper evaluates
its last stage at `t0 + h`, and `t0 + (t1 - t0) != t1` in floating
point.  `evolve_apply` now treats `h == t1 - t0` as the final step and
backs `h` off with `nextafter` until `t0 + h` can no longer pass `t1`.

The first fix alone was not sufficient: `bsimp` divides the step into
sub-steps and builds their times by repeated `t += h`, which drifts
past `t0 + h` as well.  For the final step of the reproducer the
accumulated time could dip below `t1` (for `N = 10`, to
`9.999778782798785e-13`).  `bsimp_step_local()` now computes each
sub-step time directly from `t0` and `h_total`, so the last evaluation
lands exactly on `t0 + h`.  Both the `ode-initval2` and the legacy
`ode-initval` copies carry the same code and both are fixed.

The regression test runs for every stepper: an RHS records the minimum
time seen, and the test asserts it is not below `t1` and that the final
time equals `t1` exactly.  Negative control: with the four files
reverted the new check fails for essentially every stepper.  With the
fix, `ode-initval2` is 1184/1184 and `ode-initval` is 3542/3542.

### `#39713` — secant "derivative value is not finite" — already fixed (inherited)

The report is fixed in the tree, not by this fork.  Upstream
`d50dc70e5` returns `GSL_SUCCESS` when the previous function value is
exactly zero, and `eaaae349d` replaced the divided difference
`df_new = (f_new - f)/(x_new - x)` with the algebraically equivalent
`df_new = df*((f - f_new)/f)`, which cannot divide by zero.  The posted
`gsl-roots-secant{,-v2}.patch` are against the pre-2013 line and no
longer apply; the reviewer's "skip the update for a sub-ulp step"
variant is not taken, because it would loop rather than terminate when
the step cannot move the iterate.  `ce90a6b91` added `func7` and the
vector `-pi x + e {1.5}` to `roots/test.c` as the negative control.
The reporter's own `gsl-secant64.c` (the "almost linear" case with a
`0.001*eps` offset) converges on the built DLL to
`# f(x_i) = 2.2204e-19`.  No fork change.

### `#50712` — lm+accel with a finite-difference `fvv` — rejected

The corresponding test is disabled under `#if 0` in
`multifit_nlinear/test_fdf.c` with the comment "box3d test fails on
MacOS here", and `multifit_nlinear/TODO` item 5 records it.  Enabling
it temporarily reproduces the failure on MSVC x64, so it is not
i586-only:

    FAIL: trust-region/levenberg-marquardt+accel/scale=more/
          solver=cholesky/fdfvv/box3d did not converge,
          status=exceeded max number of iterations

(the report named the svd solver; several fail).  The failing path is
the finite-difference second directional derivative
`gsl_multifit_nlinear_fdfvv()` (`h_fvv = 0.02`) with the `lmaccel`
trust region and the `box3d` problem.  The report and the upstream TODO
both leave it disabled, and no small, defensible fix was found in the
time-boxed investigation - the cause is a numerical study of the `fvv`
step and the acceleration, not a wrong-value defect.  Recorded
`[rejected]`; the `#if 0` block is left as upstream has it.

### `#30540` — convergence checks in `rk4imp`/`rk2imp` — fixed (taken later)

The patch targets the legacy `ode-initval` (v1) steppers, which run a
fixed three iterations and carry the comment "This method does not
check for convergence of the iterative solution!".  The reporter's
claim reproduces: with the Henon-Heiles program the maximum energy
error grows from `1.03e-6` (t in [0, 1000]) through `2.34e-6`
([4000, 5000]) to `4.07e-6` ([9000, 10000]) at `h = 0.1`, i.e. secular
drift rather than the bounded error a symplectic method should show.

The maintained `ode-initval2` steppers already check convergence, so
this was rejected as an enhancement to a superseded module.  It was
taken later once the filter was widened: both v1 steppers now iterate
until the stage increment satisfies
`|delta| <= GSL_DBL_EPSILON (|Y| + 1)` (cap 1000), and a step that does
not converge leaves `y` untouched and reports `yerr = GSL_DBL_MAX` so
the adaptive control rejects it.  A hard failure was rejected because
v1 `gsl_odeiv_evolve_apply()` does not retry on stepper failure and
because `gear2` uses `rk4imp` as its primer.  On the fixed-step
Henon-Heiles regression the energy error is now bounded (rk4imp about
`1.0e-9`, rk2imp about `1.6e-5`).  Description in `FORKNEWS`; commit
`fa1622e11`.

### `#30947` — fixed step size control object — fixed

A feature request for a new public object and constructor,
`gsl_odeiv_control_fixed_new()`, in the legacy v1 API.  It was taken
once the eligibility rule was widened: `ode-initval/cfxd.c` provides the
controller (`379750b93`).  The v2 interface already had
`gsl_odeiv2_evolve_apply_fixed_step()` and
`gsl_odeiv2_driver_apply_fixed_step()`.

Recorded in `FORKNEWS` under `[upstream]` (two changes) and
`[rejected]` (three reports).  Full CTest suite 56/56 on MSVC x64.


## Group P - off-track / likely reject: #42058, #45797, #59914,
## #65728

Reviewed 2026-10-04 against the built `build-cmake/gsl.dll` (MSVC x64).
This is the low-value tail: four reports that, on their face, ask for
nothing that the eligibility rule admits.  The review confirms that for
each of them, from the offline dossiers under
`temp/savannah-store/dossiers/` and the `bug-gsl` mbox threads, and
records the reasoning so they are not revisited.  **No code, test or
documentation change follows from any of them.**

The scope note: the user's working list named "43159-style items"
alongside these four.  #43159 is **not a GSL bug** - it is a GNU Octave
bug ("The Signal package (1.3.0) specgram() function is broken ..."),
closed as a duplicate of Octave #42043, and it is absent from the GSL
triage index and from every GSL source.  It is therefore dropped, not
reviewed.

Deferred for lack of offline text: none of these four is deferred;
each has a dossier.  (The four GSL bugs that were once thought to have
no local copy at all are now all covered: `50712` by Group O, and
`51104`, `53903` and `53904` by Group Q below.)

### `#42058` - GSL RSS Feed does not validate - rejected

The report is that the link labelled "Savannah GSL RSS feed" on the
GSL project page points at an ATOM feed with an XML error:

    http://savannah.gnu.org/news/atom.php?group=gsl
    line 228, column 100: XML parsing error: undefined entity

The feed is generated and hosted by **Savannah** (`news/atom.php`,
`savannah.gnu.org`); no part of it is produced by the GSL source tree.
The reporter even offers to fix it himself.  There is no GSL defect here
and nothing in this repository that could be changed: the fault, if it
persists, belongs to the Savannah `news` module or to the project's use
of an undefined HTML/XML entity in an announcement.  Savannah was
unreachable during this review (`news/atom.php` returned a transport
error), which is consistent with the meta/service nature of the report.
Rejected - not a GSL defect; outside the repository.

### `#45797` - "Possible problem with LAPACK Fortran routine ZGESVD" - rejected

The reporter built a 4x4 complex SVD test with **LAPACK's `ZGESVD`** and
saw the imaginary part of one left singular vector differ between two
Ubuntu kernels (`2.6.31-16` "correct" and `2.6.32-40` "incorrect").  The
report is explicitly about the LAPACK/ZGESVD behaviour bundled by the
distribution, not about a GSL routine.

GSL does not expose a complex SVD at all.  The public interface has only
the real factorisations - `gsl_linalg_SV_decomp`, `SV_decomp_mod`,
`SV_decomp_jacobi` in `linalg/gsl_linalg.h`; there is no
`gsl_linalg_complex_SV*` and no `zgesvd`/`ZGESVD` reference anywhere in
the sources (checked by search).  The reporter's program links LAPACK
directly; there is no GSL code path to reproduce, and a change of LAPACK
or kernel behaviour is not something GSL controls.

Rejected - external LAPACK/distribution report; no GSL code involved.
No reproduction is possible because the named routine does not exist
here.  (A complex SVD would be a feature request, and one is in fact
already filed separately as #66574/#66575.)

### `#59914` - "Native build of GSL-2.5 on windows 10" - rejected

The report is a *success story*, not a defect: the submitter cross-built
GSL 2.5 DLLs on Ubuntu with MXE and offers to document the recipe on a
git page, asking whether it would be worth linking from the GSL website.
It carries no patch and reports no wrong result.

The concrete need - a maintained native Windows build - is already met
in this fork by the CMake build, and by the Windows CI and DLL workflow:
`CMake.md` documents the MSVC build ("Windows with MSVC, tested with
Visual Studio 2022 Build Tools, x64"), `.github/workflows/cmake.yml`
runs `windows-msvc` shared and static, and
`.github/workflows/windows-dll.yml` packages a Windows x64 build.  The
remaining part - a link from the upstream project page to a
community write-up - is a website/documentation matter outside this
repository and not a bug fix or doc *correction* in the tree.

Rejected - no defect and no actionable change; the Windows building path
exists in the fork already.

### `#60026` - "Incorporate MIXMAX random number extension into GSL" - fixed

A request to add the MIXMAX family of PRNGs to GSL, citing its adoption
in CLHEP/Geant4, CMS simulation and a NASA neutrino-telescope study.  The
2025 follow-up locates a reference implementation in C++ (CLHEP
`Random/src/MixMaxRng.cc`, GPLv3) and notes it would need translating to
C; a further reply asks for TestU01/PractRand evidence and states none of
the MIXMAX variants are in the reviewer's own Dieharder quality table.

Taken once the rule was widened (like the neighbouring PCG request,
#55965).  The N = 17, one-parameter variant over GF(2^61 - 1), the
configuration CLHEP selects by default, is added as `gsl_rng_mixmax17`
(260557ee9).  The native draw is 61 bits; `gsl_rng_get` returns the low
32 bits, so the value fits `unsigned long` on LLP64 and LP64, while
`gsl_rng_uniform` uses the full 61-bit value.  Seeding uses the reference
`seed_spbox` procedure, which needs no 128-bit multiply, so the
portability concern that kept the ticket in the backlog does not arise.

Verified against ROOT's independent bundled MIXMAX implementation with
N = 17 and the `seed_spbox` path (embedded vectors for seeds 1, 42 and
12345 after 1, 10^4 and 10^6 draws), with a negative control on
SPECIALMUL.  MIXMAX is linear over GF(2^61 - 1) and so fails
modular-rank tests, as MT19937 and PCG do; TestU01/PractRand are not in
the suite.  Full CTest suite 56/56 on MSVC x64.

### `#65728` - "Add sparse functionalities" - partly fixed

A feature request forwarded by the maintainer, asking for complex sparse
support: scaling a complex sparse matrix by a *real* vector, complex
sparse matrix-vector multiplication, in-place unpacking of a complex
vector, resize operations on blocks/vectors, and an accumulating variant
of `_set()` that does `*ptr += x` rather than replacing.

Moderate parts of this already existed: the complex sparse type has
`gsl_spmatrix_complex_scale`, `_scale_columns`, `_scale_rows`, `_memcpy`,
`_dense_add`/`_dense_sub` and `add_to_dense`, and the complex block/vector
types have `alloc`/`calloc`.  Two of the missing pieces were taken once
the eligibility rule was widened (`f218d8c83`): the complex sparse
matrix-vector product `gsl_spblas_zgemv()` and the accumulating
`gsl_spmatrix_*_set_add()`.  The `resize`, real-vector scaling and
in-place unpack parts remain in the backlog.

Recorded in `FORKNEWS` under `[rejected]` (one combined entry).
No code change, so the full CTest suite is unchanged at 56/56 on MSVC
x64.


## Group Q - the last open items: #51104, #53903, #53904

Reviewed 2026-10-04 against the built `build-cmake/gsl.dll` (MSVC x64).
These three were the last open reports without a recorded verdict; each
has an offline dossier under `temp/savannah-store/dossiers/`, so the
lack of a local copy that had deferred them no longer applies.  All
three are rejected; **no code, test or documentation change follows**.
Recorded in `FORKNEWS` as one combined `[rejected]` entry.

### `#51104` - gsl_permutation_next efficiency - rejected

The report (2017, Performance) offers a more literal translation of
Knuth's Algorithm L (AoCP volume 4, page 319) and reports a speed-up,
linking a gist; the single follow-up adds that the second comparison in
the selection loop, `(p->data[j] > p->data[i]) && (p->data[j] <
p->data[k])`, is redundant because the loop that finds `i` already
establishes `p->data[j] < p->data[k]`.

No defect is alleged and the current implementation is correct;
`permutation/test.c` walks the full lexicographic order and passes.  The
proposal is a performance-only change, which the eligibility rule
excludes outright; the cache-miss loop reorder #54925 was rejected on
the same ground.  Rejected - performance-only, out of scope.

**Superseded (commit `c382ed86e`).**  The selection loop was changed to
the single backward scan the report proposes; the generated permutations
are unchanged.  `#21833` is the same function and is closed with it.
Applied after the filter was widened; see Group V at the end of this
file.

### `#53903` - Test failure with gsl_sf_synchrotron_1_e on x86 - rejected

The report (2018, Runtime error) is that `gsl_sf_synchrotron_1_e(0.01)`
fails its test on an i386 build configured for Pentium MMX
(`-march=pentium-mmx -mtune=pentium-m`, musl libc 1.1.19, gcc 6.4):

    expected: 4.4497250411421063e-01
    obtained: -1.8137993642342178e-02 +/- 1.2082330897339797e-17   (rel=6.66e-16)

The reported numbers pin the fault to the toolchain rather than to GSL.
`-1.8137993642342178e-02` is exactly `-c0*x` with `c0 = M_PI/M_SQRT3`
and `x = 0.01`, and the reported error `1.2082330897339797e-17` is
exactly `3*eps*c0*x`.  In `synchrotron.c` the value for `x = 0.01` comes
from the `x <= 4` branch,

    result->val  = px * c1.val - px11 * c2.val - c0 * x;

with `px = pow(x, 1.0/3.0)` and `px11 = px^11`.  The reported value and
error are what this branch produces only when `px` is zero, which drops
both Chebyshev terms and leaves the error as
`c0*x*eps + 2*eps*|val| = 3*eps*c0*x`.  The Chebyshev coefficients
themselves are correct: evaluated independently (Python, with the exact
coefficients from `synchrotron1_cs`/`synchrotron2_cs`) the branch gives
`0.4449725041...`, matching the test's expected value, and the specfunc
suite passes on MSVC x64 (`specfunc_test.exe`, exit 0).

So a 32-bit x87 target's libm `pow` returned zero for
`pow(0.01, 1.0/3.0)`; the same class of toolchain fault as the MinGW
`sin`/`cos` argument reduction already recorded under Group H.  There is
no 32-bit x87 host here to reproduce on, and correct code is not changed
for an unreproducible libm fault.  Rejected - toolchain, not a GSL bug.

### `#53904` - Bug gsl_matrix_complex_set - rejected

The report (2018) is that passing `gsl_complex_rect(1., 0.)` directly to
`gsl_matrix_complex_set(m, 1, 1, ...)` corrupts the matrix, while
assigning it to a `gsl_complex` variable first works.  The program
allocates the 4x4 complex matrix with `gsl_matrix_complex_alloc` (which
does not zero it), sets only element `[1][1]` - inside a 16-iteration
loop - and runs `gsl_eigen_herm` on it, printing neither the matrix nor
the eigenvalues; the "wrong data array" it prints is the program's own
input array.

The two forms are equivalent C.  `gsl_complex_rect` is an `INLINE_FUN`
returning a `gsl_complex` by value, and `gsl_matrix_complex_set` takes
its last argument by value as well, so storing the result in a
`gsl_complex` first changes nothing.  The difference the reporter saw is
the uninitialized remainder of the matrix making the eigenvalues
garbage.  Rejected - not a bug.

Savannah bugs #51104, #53903 and #53904.


## Group R - the last bug-shaped reports: #42830, #68068, #47348,
## #68398, #68611 (fixed) and #40116, #45099, #68663, #68704 (rejected)

Reviewed 2026-10-04 against the built `build-cmake` library.  Nine
reports; five are genuine defects and are fixed, four need no code
change.

### Fixed

* **#42830 - bspline breakpoint validation.**  `gsl_bspline_init_augment`
  wrote the caller's `tau` straight into the knots with no monotonicity
  check, so a decreasing breakpoint vector gave a non-monotone knot
  sequence and meaningless (negative) basis functions.  The guard is
  non-decreasing rather than strictly increasing, because
  `gsl_bspline_init_greville` legitimately produces repeated interior
  breakpoints; a few-ulp reversal from its least-squares solve is
  tolerated (observed `4.0000000000000018, 3.9999999999999978` at
  `k=4, nbreak=5`).  Commit `74cace84f`.  Negative control: the new
  `bspline/test.c` vector reports `GSL_SUCCESS` with the guard removed.

* **#68068 - quad_golden `*f_upper` typo.**  In the branch that moves
  `x_lower` up to `x_m`, the routine stored `f_m` in `*f_upper` instead
  of `*f_lower`.  Iteration was unaffected (it reads only the `x`
  bounds) but `gsl_min_fminimizer_f_lower` returned 0.65 where
  `f(x_lower) = -0.116` on the `func4` vector.  Commit `4260778fe`.
  Negative control: an added per-iteration invariant check in
  `min/test.c` fails with dozens of mismatches without the fix.

* **#47348 - floor(x+0.5).**  Swept all 31 live sites to `rint`.  The
  `hyperg_1F1.c` "b-a is an integer" test was additionally missing its
  `fabs()`, so it accepted any `b-a` with fractional part > 0.5; a
  14950-point sweep showed the affected range takes the same
  `hyperg_1F1_U` fallback either way, so no value changed.  Two sites
  are deliberately *not* rounded to nearest:
  `specfunc/coulomb.c`'s `lam_F` reduction is a round-half-up (keeps
  `lam_min` on the well-conditioned `-1/2` branch; `rint(0.5)=0` would
  leave `lam_min=+1/2` and return `GSL_EDIVERGE`), left as-is with a
  comment and pinned by a new vector.  Commit `b681165e1`.  Observed
  effect: `gsl_sf_bessel_Ynu_e(0.5, 1.0)` moved by less than an ulp
  (both accurate), added as a vector; no other tested value changed.

* **#68398 - eta_int leading term.**  `2^(1-n)` is exactly
  representable, so `gsl_ldexp` beats `gsl_sf_exp_e`: at `n = -101` the
  value is correctly rounded (2.75e-16 vs a 40-digit reference, from
  5.6e-14) and the reported error halves.  Commit `f774247e7`.
  Negative control: the new error-bound check reads 5.56e-14 with the
  exponential restored.

* **#68611 - all-empty histogram pdf.**  Zero mean made every cumulative
  sum `NaN`, so `gsl_histogram_pdf_sample` returned `NaN` silently.
  Now `GSL_EDOM`, matching the existing negative-bin check.  Commit
  `b60a47106`.  Negative control: the trap test fails without the guard.

### Rejected

* **#40116 - "possible error in integration routines".**  Faithful
  QUADPACK restructuring: `large_interval` is `if(|b-a|>small) go to
  90` with `level[i] < maximum_level` as the exact equivalent, and
  `increase_nrmax()` returning 1 is the `go to 90` (so the loop does
  cycle).  No buggy case.  Not a bug.

* **#45099 - "wrong BFGS update".**  Expanding the dense update with
  `H0 = I` gives `Hg = g - A s - B y`, `B = s.g/s.y`,
  `A = -(1+|y|^2/s.y)B + y.g/s.y` - exactly the code.  The Hessian is
  implicit in `p`, not stored.  Not a bug.

* **#68663 - remove `GAMMA_INC_A_0`.**  One-line alias for
  `gsl_sf_expint_E1_e`; pure refactor, out of scope.

* **#68704 - `gsl_histogram_variance`.**  New public API; the proposed
  rename of `gsl_histogram_sigma` is a breaking change.  Out of scope.

Method notes:

* Two negative controls were checked per new vector: reverting the
  source fix (the vector fails) and reverting the test (the suite
  passes), to be sure a failure came from the change under test and not
  from a neighbouring test-vector edit.
* A new Coulomb vector briefly corrupted the existing
  `gsl_sf_coulomb_wave_FG_e(-50, 1000, 0, 0)` test by leaving `lam_F`
  set from a preceding case; the fix restored the `lam_F = 0` and `k_G
  = 0` assignments the block had relied on, and the full `specfunc`
  suite then reported `59621340/59621340`.
* Full `ctest` after all changes: 56/56 on MSVC x64.

Savannah bugs #40116, #42830, #45099, #47348, #68068, #68398, #68611,
#68663 and #68704.


## Group U - the remaining feature/perf/API requests: closed as rejected

Reviewed 2026-10-04 from the offline dossiers under
`temp/savannah-store/dossiers/` and the `bug-gsl` mbox threads.  This is
the tail the eligibility rule does not admit, recorded in one place so
the index has no un-triaged item left.  **No code, test or documentation
change follows from any of them** at the time of writing; several were
taken later once the filter was widened (the six special cases, `#68098`,
`#24871` and some performance items; see `FORKNEWS`).  The detail is in
`FORKNEWS` under "meta: the remaining feature, performance and API
requests"; the groups below are the shape of the request, not a per-bug
essay.

* New distributions/generators (new public API): `#24252`,
  `#59900`, `#66767`, `#66775`, `#66816`, `#66949`.
  **Partly reversed:** `#24871` (complex exponential integrals E_n) was
  taken later as a port of Amos Algorithm 683 (`c01e6b922`); see
  `FORKNEWS`.
* New algorithms/solvers: **all taken** — `#66573`
  (`gsl_vector_complex_conjugate()`), `#68367` (inverse Jacobi elliptic
  functions), `#60457`/`#65932` (complex tridiagonal solvers), `#66574`
  (complex Householder right, Givens, bidiagonal decomposition and SVD),
  `#66575` (complex SVD solve and modified decomposition), `#66695`
  (Feagin/Verner ODE steppers), `#57173` (complex zeta, eta and Hurwitz
  zeta) and `#32257` (Gauss-Lobatto/Radau/Clenshaw-Curtis/Fejer
  quadrature) were taken once the filter admitted new public API; see
  `FORKNEWS`.
* New special cases or API extensions: `#41527`, `#66800`,
  `#66834`, `#66842`, `#66850`, `#66922`, `#67359`, `#67774`, `#68098`.
  **Partly reversed:** the six among them (`#66800`, `#66834`, `#66842`,
  `#66850`, `#66922`, `#67359`) were taken as correctness fixes,
  `#68098` was taken later as an argument-range fix (`6da0c08f3`) and
  `#45782` (configurable finite-difference Jacobian step) was taken as
  new public API (`8abbe2c30`).  `#41527` was taken later as a bug fix
  once the filter admitted contract changes to existing API: the simplex
  minimizers propagate `GSL_EBADFUNC` instead of masking a non-finite
  objective as `GSL_EFAILED` (`f007aa42c`); it is no longer an API
  extension.  `#67774` remains out of scope.
* Breaking/behaviour change: `#68549`.
* Documentation/test/refactor with no patch, or a reorganisation:
  `#66742`, `#66826`, `#66844`, `#66874`, `#66877`, `#66880`, `#66886`.
* Performance-only (changes the floating-point result, no wrong value
  alleged): `#21828`, `#21833`, `#31109`, `#40092`, `#51104`.
  **Partly reversed:** `#21833`, `#40092` and `#51104` were taken later
  as result-preserving fixes (`c382ed86e`, `0697ba772`); `#54925` was
  taken too (`39cde03e2`).  `#21828` and `#31109` remain out of scope.
  See Group V at the end of this file.
* Not a defect: `#47402` (Mathieu design discussion).

Two of these deserve a note:

* `#66874` is **superseded**: the binary search tree module is already
  documented in the fork (`20eda0316`), so the posted `doc_bst.rst` is
  not needed.
* `#66742` (GAMS classification) is a 177-file documentation sweep, not
  a correction of anything wrong, so it is a feature, not a doc fix.

The four reports once listed here as still under review - `#21831`
(Levy skew for alpha < 1), `#25320` (Fresnel extension), `#29834` (BLAS
wrapper argument checking) and `#34361` (`bspline_knots_greville`
constrained least squares) - are resolved in Group C below.

Savannah bugs #21828, #21833, #24252, #24871, #31109, #32257, #40092,
#41527, #45782, #47402, #51104, #57173, #59900, #66573, #66695, #66742,
#66767, #66775, #66800, #66816, #66826, #66834, #66842, #66844, #66850,
#66874, #66877, #66880, #66886, #66922, #66949, #67359, #67774, #68098,
#68367 and #68549.


## Group C - the last numerical reports: #21831, #25320, #29834, #34361

Reviewed 2026-10-04 against the built `build-cmake/gsl.dll` (MSVC x64)
and the offline dossiers under `temp/savannah-store/dossiers/`.  These
were the four reports still open at the end of Group U; three describe
no defect in the current tree (or were fixed long ago) and one is a real
hazard that the `init_augment()` guard already rejects.

### `#21831` - Levy random number generator for alpha < 1 - not a bug

The report (2007) has two claims.  First, that `gsl_ran_levy_skew()`
with `beta = 0` does not produce a Levy variate but a "linear" density.
The symmetric case has delegated to `gsl_ran_levy()` since the skew
function was introduced (`8262207e8`, 2001); driving both with the same
seed gives bit-identical streams.  Second, that `gsl_ran_levy()` loses
accuracy for `alpha < 1`, requiring a sum of many variates, and is worse
still below `0.3`.

Both claims were checked against the characteristic function rather than
against another sampler.  The transform is the Chambers-Mallows-Stuck
one, whose exact characteristic function is in the source comment; for
the symmetric case,

    E[exp(i t X)] = exp(-|c t|^alpha)

and for the skew case,

    E[exp(i t X)] = exp(-|t|^alpha (1 - i beta sign(t) tan(pi alpha/2)))
                    (alpha != 1)

Calling the built DLL through ctypes and averaging `cos(tX)`, `sin(tX)`
over 400000 draws gives agreement to the Monte-Carlo floor
(`~2e-3`, i.e. `1/sqrt(N)`):

    alpha=0.5, beta=0:  E[cos(tX)] vs exp(-|t|^0.5)   diff <= 2.1e-3
    alpha=0.8, beta=0:  E[cos(tX)] vs exp(-|t|^0.8)   diff <= 2.1e-3
    (alpha,beta)=(0.5,0.5):  matches both components  diff <= 3e-4
    (alpha,beta)=(1.0,0.5):  matches both components  diff <= 3e-4
    (alpha,beta)=(1.3,0.7):  matches both components  diff <= 3e-4

The reporter compared a `10^6`-sample histogram against a numerical
integration of the oscillatory characteristic function; that reference
is the inaccurate part for heavy tails.  Rejected - not a bug.  The
test suite does not exercise `alpha < 1`, but a statistical regression
test would be seed- and tolerance-sensitive, so none is added.

### `#25320` - Import fresnel, bugs on GSL Extension Fresnel - rejected

The report is a maintainer note to import Andrew Steiner's Fresnel
extension, carrying Toshiro Ohsaki's observation that the extension
returns the wrong sign for negative `x`.  The observation is right, but
there is no Fresnel code in GSL to fix: a repository-wide search finds
`fresnel` only in `TODO` and `specfunc/TODO`, both asking for the
integrals to be added.  The extension lives outside the tree.  Adding
Fresnel integrals (and fixing their negative-`x` branch) is a new
special function - new API and a new algorithm - which the eligibility
rule excludes.  Rejected / deferred as a feature.

### `#29834` - insufficient argument checking in blas wrapper - already fixed

The report (2010) is the maintainer's note that the CBLAS routines do
not validate their arguments and that the checking could be shared
through macros.  The patch attached to the report
(`error_cblas_v2.h`) was integrated shortly afterwards: `cblas/error_cblas.h`
holds the `CHECK_*` primitives and `error_cblas_l2.h` /
`error_cblas_l3.h` the per-routine `CBLAS_ERROR_*` macros, called from
the `source_*.h` kernels (commits `3b3c12ef1` "added error checking",
`36434582b` "revised cblas error checking patch from jgpallero",
`f30cb5cd6`, `e9a1275f0`).  The `gsl_blas_*` wrappers already perform
their dimension checks and return `GSL_EBADLEN` (for example
`gsl_blas_dgemv()` in `blas/blas.c`).  The remaining wish - unify the
two layers through the macros - has no failing case attached and is an
open-ended refactor.  Rejected / already addressed; no fork change.

### `#34361` - gsl_bspline_knots_greville needs constrained least squares - hazard fixed

The function selects its breakpoints by an unconstrained linear
least-squares solve for the target Greville abscissae, without enforcing
that they be non-decreasing.  Replicating `greville.c`'s linear algebra
in NumPy, 86448 of about 180000 random monotone requests (any
`k = 4..8`, any `nbreak = 4..12`) produced a non-monotone breakpoint
vector, e.g. `k = 8, nbreak = 4` with abscissae
`{0.0745, 0.441, 0.462, 0.5876, 0.6113, 0.6795, 0.6804, 0.7269, 0.7329,
0.9872}` gave `{0.0745, 2.5803, -0.8523, 0.9872}`.

The guard added to `gsl_bspline_init_augment()` in `74cace84f`
(Savannah bug #42830) already rejects this case with `GSL_EDOM`, before
the workspace is written, and its comment names this caller.  The
missing pieces were a Greville-specific test and documentation, both
added in `c8543330d`:

  * `bspline/test_greville.c`: a `k = 4`, `nbreak = 5` call with the
    abscissae `{0, sqrt(1/6), ..., 1}` (unconstrained breakpoints
    `{0, 1.0611, 0.6091, 0.6767, 1}`) must return `GSL_EDOM`.  Negative
    control: with the `init_augment()` guard disabled, `status = 0`.
  * `doc/bspline.rst` + `doc_texinfo/bspline.texi`: the
    `gsl_bspline_init_greville()` / `gsl_bspline_knots_greville()` entry,
    commented out since the bug was filed, is restored with the
    `GSL_EDOM` exception and the workspace-unchanged guarantee.

The full LSI solver the report proposes (Lawson-Hanson NNLS/LDP) is a
new algorithm and remains out of scope; the source still marks the
routine "Limited function".

Full `ctest` 56/56 on MSVC x64.

Savannah bugs #21831, #25320, #29834 and #34361.


## Group V - the performance-only items, re-examined

Reviewed 2026-10-04 against the built `build-cmake` library (MSVC x64).
The eligibility rule had set every performance-only report aside.  The
filter was widened for changes that are either numerically neutral
(bit-for-bit identical results) or remove provably redundant work, and
that need no interface change.  Four of the six backlog items qualify
and are applied; the other two do not.  No platform-conditional
tolerance or benchmark result is used to justify a change.

### `#54925` - source_gemm_r loop reordering - applied

The posted patch swaps the `i`/`k` loop nesting in the two `NoTrans`
branches of `cblas/source_gemm_r.h` and the `i`/`j` nesting in the two
`Trans` branches, so a row of `C` is traversed once and stays resident
while `F` and `G` are read.  The earlier rejection rested on "reordering
the accumulation changes the floating point result", which is wrong:
for each `C(i,j)` the terms are still added in ascending `k` order, and
only the order in which the independent `(i,j)` pairs are visited
changes.  The result is bit-for-bit identical, so the change is
numerically neutral and was taken.

Verification: the existing 104 `cblas` gemm cases already span both
storage orders and all four transpositions; a new triple-loop reference
test, `cblas/test_gemm_loops.c`, checks every combination plus non-`1`
`alpha`/`beta` and passes.  Full `ctest` 56/56.  Commit `39cde03e2`.

### `#51104` + `#21833` - gsl_permutation_next - applied

The selection loop walked the whole monotone suffix and kept the
smallest candidate, re-testing it with a comparison already implied by
the loop that located `i`.  Both `next` and `prev` now scan from the
right and stop at the first qualifying element, i.e. Knuth's Algorithm
L step L3.  The suffix is monotone, so the same `k` is chosen and the
output is unchanged.

`#21833` (2007) is the same function.  Its "quadratic worst case" claim
does not hold: summed over all `n!` permutations the reversal work is
amortised `O(1)` per call.  The regression test was raised from the
existing 5-element vectors to an exhaustive `8! = 40320` walk in both
directions against an independent reference, which also serves as the
record for `#21833`.  Full `ctest` 56/56.  Commit `c382ed86e`.

### `#40092` - false-position excess evaluations - applied

`falsepos_iterate()` evaluated `f` at the linear interpolation point and
then, when the step failed to halve the bracket, again at the bisection
point.  When the interpolation collapses onto an endpoint (an almost
linear function with an endpoint within rounding of the root) the first
evaluation repeats a known value, so every iteration costs two
evaluations instead of one.

The routine now reuses `f_lower`/`f_upper` for `x_linear == x_left` /
`x_right` and evaluates only otherwise.  This is the reporter's "quick
fix", not the complete fix (which needs the tolerance passed into the
iterate routine and an interface change) - but it is exact, not a
tolerance heuristic: the value skipped is the one the endpoint already
carries.

Verification: the report's 64-bit test function drops from 98 to 52
evaluations; `roots/test.c` counts evaluations on the same function and
requires at most 60.  Negative control: with the skip removed the count
is 98 and the test fails.  Full `ctest` 56/56.  Commit `0697ba772`.

### `#21828` (lmsder) and `#31109` (bsimp) - still out of scope

Neither admits a result-preserving change from the current source:

* `#21828` is a 2007 comparison against netlib MINPACK with no patch and
  no wrong value.  The likely costs are the column-oriented Householder
  work in `linalg/qrpt.c`/`householder.c` on GSL's row-major storage and
  a duplicated `compute_gradient_direction()` in `lmiterate.c`/
  `lmpar.c`; confirming either needs a profiling build and an
  independent reference, which was not done here.
* `#31109` needs the requested tolerance inside `bsimp_apply()`, which
  is an interface change; there is no internal signal to key the order
  on.  Left as a design item.

Savannah bugs #54925, #51104 and #21833, and #40092.


## Open defect: `gsl_sf_hyperg_1F1_int_e` for negative integer b

Reported through Savannah bug #66826 but not fixed there (the change for
that bug is test-only; see `FORKNEWS`).  Recorded here with the full
analysis so the next reader does not have to re-derive it.

For a nonpositive integer `a` and a negative integer `b` with `b <= a`,
`1F1(a,b,x)` is the terminating polynomial of DLMF 13.2.2,
`sum_{k=0}^{-a} (a)_k/(b)_k x^k/k!`.  The library disagrees:

* `gsl_sf_hyperg_1F1_int_e` returns `exp(x)` when `a == b` (line 1810);
  for `a` a negative integer the true value is the truncated
  exponential, not `e^x`.
* otherwise `hyperg_1F1_ab_negint` maps `x < 0` to `x > 0` with the
  Kummer transformation `exp(x) M(b-a,b,-x)`.  That map is not the
  terminating polynomial and disagrees with it by an amount that grows
  without bound with `|x|`; e.g. at `(-10,-20,-100)` the map gives
  `1.64e-35` against a true `4.93e7`.

### Why the obvious fixes fail

* **Direct terminating polynomial (Horner).**  Correct in principle but
  catastrophically ill-conditioned for large degree and `|x|`.  In the
  same double recurrence the library uses:

      M(-100,-200,-100) = 58.6      (true 4.34e-25)
      M(-1000,-2000,-100) = 1.93e4  (true 1.03e-22)
      M(-10,-100,-100)    rel 6.9e-9

  This is genuine cancellation, not a coding error; some values cannot
  be recovered in double by this route at all.

* **Kummer to positive argument.**  Not an identity for the terminating
  branch.  Agreement is coincidental and parameter-specific:

      (-10,-20,-1)     rel diff 1e-25
      (-10,-20,-10)    rel diff 3.4e-4
      (-10,-20,-100)   rel diff 1.0
      (-100,-200,-100) rel diff 1e-29
      (-10,-100,-100)  rel diff 0.98
      (-3,-5,-1)       rel diff 1.7e-4

  There is no simple validity test to key a fallback on.

### Consequence for the existing tests

Three vectors in `specfunc/test_hyperg.c` assert the Kummer value, which
equals the terminating series for some parameters and not others:

    (-10,-20,-10)   0.00357079636732993491  vs true 0.00357201182794578047
    (-10,-20,-100)  1.64284868563391159e-35 vs true 49303272.262405369733
    (-10,-100,-100) 8.19512187960476424e-09 vs true 4.4148250205403363517e-7

They pass today because the library returns the same wrong number, and
must be corrected when the library is fixed.

A correct fix needs a stable evaluation for these polynomials (a scaled
recurrence, a Laguerre/Whittaker route that handles the negative integer
parameters, or an error-bounded series), plus re-derived expected values
and an independent cross-check over a grid.  That is a separate task.


## Group W - two already-classified items taken (2026-10-06)

Reviewed against the built `build-cmake/gsl.dll` (MSVC x64 Release).
Both were earlier set aside (`#67774` as an API extension, `#21828` as
performance-only) and are now taken: `#67774` turned out to be a defect
in behaviour the library already claimed, and `#21828` admits a
result-preserving change.

### `#67774` - the arctangent integral's error bar for x < 0 - fixed

The 2025 report frames this as extending the domain, but the domain was
never restricted: the manual defines
`AtanInt(x) = int_0^x arctan(t)/t` with an empty `Domain`, and the code
already took negative arguments.  What was wrong was the error bar.  In
the `|x| <= 1` branch the series value and error come from a Chebyshev
fit in `x*x`,

    result->val = x * c.val;
    result->err = x * c.err;

so the error is multiplied by the signed `x` (`ax` was already computed
for the branch tests).  Measured with the built DLL through ctypes:

    x = -1    val = -0.91596559417721901   err = -2.93e-17
    x = -0.5  val = -0.48722235829452237   err = -6.44e-18
    x = -0.1  val = -0.099889286860336185  err = -1.39e-18
    x = -2    val = -1.5760154034463234    err =  8.15e-16   (branch |x| > 1)

The value is correct (the function is odd); only the error's sign is
wrong.  The fix is `err = ax * c.err`; the `|x| > 1` branch already used
a non-negative expression.

Verification: four vectors (`-0.1`, `-1.0`, `-2.0`, `0.0`) added to
`specfunc/test_sf.c` at `TEST_TOL0`, mirroring the existing positive
values (which the report's `-1`/`-2` values also match).  Negative
control: with the one-line change reverted, `specfunc_test` fails on
`(-0.1)` and `(-1.0)` with both "reported error negative" and "value not
consistent within reported error"; `(-2.0)` is unaffected because it
takes the `|x| > 1` branch.  The expected values were checked
independently against mpmath `quad(atan(t)/t, [0, x])` at 40 digits.
Full ctest 56/56.  Commit `4220f05d9`.

### `#21828` - the lmsder Householder cost - result-preserving rewrite

The 2007 report compares `gsl_fdfsolver_lmsder` with netlib MINPACK and
names the Householder transform (about 50 % of the run time); it has no
patch and no wrong value.  Profiling the current source confirmed the
location but also found a second, larger cause.

`gsl_linalg_householder_hm()` (`linalg/householder.c`) applied
`H = I - tau v v^T` to the target matrix **column by column**.  GSL
matrices are row-major, so for a fixed column the inner loop over rows
stepped with the row stride.  On the CMake/MSVC build the loop was worse
still: `HAVE_INLINE` is undefined (`build-cmake/config.h`), so
`gsl_matrix_get()`/`set()` and `gsl_vector_get()` were real out-of-line
function calls in the innermost loop.  `householder_hm` is the dominant
cost of `gsl_linalg_QRPT_decomp()`, which lmsder re-enters on every
successful iteration (`multifit/lmset.c`, `multifit/lmiterate.c`).

The non-BLAS branch now walks a block of 64 columns at a time using
direct row pointers.  The dot product still accumulates in increasing
row order and each element is updated with the same expression
(`Aij - tau*vi*wj`), so the output is bit-for-bit identical; only the
order in which independent columns are visited changes.  The `USE_BLAS`
branch is untouched (the fork's builds never define `USE_BLAS`).

Measured on a single `gsl_linalg_QRPT_decomp` of a 550 x 25 matrix
(MSVC x64 Release, best of 5 x 400 calls):

    original (column order, accessors)  1.43 ms
    blocked row order, accessors        0.96 ms
    blocked row order, direct pointers  0.21 ms

Bit-identical check: the harness drives
`gsl_linalg_QRPT_decomp` over six shapes (`5x4`, `20x10`, `37x36`,
`64x30`, `100x3`, `550x25`) and three seeds, and compares the raw bytes
of the full matrix, `tau`, the permutation and the column-norm vector
between the pre-change and post-change DLLs; they are identical.  Full
ctest 56/56.  A numerically neutral kernel has no negative-control
vector to add.  Commit `cd9e31326`.

Savannah bugs #67774 and #21828.


## Group W - the new distributions: five taken, the gamma tail re-rejected

Reviewed 2026-10-06 against the built `build-cmake/gsl.dll` (MSVC x64).
The eligibility filter was widened to admit new public API, so the
"new distributions/generators" group from Group U was revisited.  Five
were implemented and one was re-rejected on correctness grounds.

### `#66949` - Erlang cumulative distribution - taken

`randist` had `gsl_ran_erlang()` but `cdf` had no Erlang cdf.  Added
`gsl_cdf_erlang_P/Q(x, k, lambda)` as the gamma cdf with shape `k` and
scale `1/lambda`.  Vectors are the exact Poisson sum
`1 - sum_{i=0}^{k-1} exp(-lambda x)(lambda x)^i/i!`; reversing the
delegation fails them.  Commit `3d607fce8`.

### `#66816` - Nakagami distribution - taken

Added `gsl_ran_nakagami()`/`_pdf()`; a Nakagami variate is the square
root of a `Gamma(mu, omega/mu)` variate.  Pdf vectors were checked
against scipy/mpmath, and the sampling against the pdf and the analytic
tail.  The consuming tests were placed after every other distribution
test: an earlier placement shifted the shared RNG stream and made the
multivariate Gaussian test fail at `p = 0.027`, which is the same
stream-sensitivity the `beta_small` comment already warns about.
Commit `999d53bbe`.

### `#59900` - truncated normal distribution - taken, self-contained

The posted patch called `gsl_cdf_*` from `randist/gauss.c`.  That would
make `randist` depend on `cdf` while `cdf` already depends on `randist`
(`betainv.c`, `hypergeometric.c`, `tdistinv.c`), and the patch fixed
only three of the affected autotools test link lines (not `cdf/test`
itself).  The implementation here keeps `randist` self-contained: the
generator rejects from the untruncated Gaussian and the pdf normalizes
with `gsl_sf_erfc`.  The cdf side reuses `gsl_cdf_gaussian_*` inside the
`cdf` module.  The Burkardt P/Q/Pinv/Qinv table was reproduced
independently with mpmath; the P/Q vectors use `TEST_TOL2` because the
difference of two Gaussians loses a few ulp.  Commit `76f3b1db8`.

### `#66767`, `#66775` - inverse binomial and Poisson - taken

Added the discrete quantiles as `double`, smallest `k` with
`P(k) >= P` (or `Q(k) <= Q`), by bisection (the Poisson range is
doubled first).  The report's Poisson prototype was a copy of the
binomial one; the real signature is `(P, mu)`.  Vectors are R's
`qbinom`/`qpois` plus exact small-case sums; shifting the result by one
fails them.  Commit `e88a81c75`.

### `#24252` - gamma tail distribution - rejected (again, for a new reason)

The 2008 attachment `gamma_tail_jpl_080908.c` was recovered from the
Wayback Machine.  It is a left-truncated gamma sampler
`devroye_tail_gamma_ran(rng, tail, a)` with unit scale, proposed again
in 2026 as a replacement for the `a < 1` branch of `gsl_ran_gamma`.

It is numerically wrong.  For `a < 1` the proposal is `Y = tail +
Exp(1)` and the target is proportional to `Y^(a-1) e^(-Y)`, so the
acceptance probability is `(tail/Y)^(1-a)`; the code uses
`(a/Y)^(1-a)`, i.e. the constant `a` where it must be `tail`.  The two
agree only when `tail = a`.  Simulated over 200000 draws:

    a = 0.5, tail = 0: sample mean 0.7476   Gamma(0.5) mean 0.5000
    a = 0.8, tail = 0: sample mean 0.9177   Gamma(0.8) mean 0.8000

(the `a >= 1` branches do match).  The `tail = 0` case the 2026 comment
proposes is exactly where the branch is meaningless, so the claimed
speed-up cannot hold.  No change; a correct truncated-gamma generator
would be separate work.  Recorded as `[rejected]` in `FORKNEWS`.

Savannah bugs #66949, #66816, #59900, #66767, #66775 and #24252.

