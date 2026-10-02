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

Adds four `test_hstaxe_1..4` functions to `interpolation/test.c`.  Two of
the four are near-duplicates of the others, all four leak the `gsl_interp`
and `gsl_interp_accel` they allocate, two declare unused variables
(`test_table`, `min_size`, `wav_min`, `wav_max`), and `test_hstaxe_2`
asserts the result is 0.0.  Low value as posted; if the coverage is
wanted, one clean test is better than four.

### `#64851` — gsl-config exit status

Superseded.  The one-line change (`usage; exit 1` → `usage 1`) is already
equivalent to what this fork's own `gsl-config` work does with
`usage 1` plus `exit $1`.


## Not yet reviewed

Everything on the curated list has now been reviewed; what remains are
feature-shaped patches excluded by choice (see the eligibility rule in the
discussion): `#66573`, `#66574`, `#66575` (complex SVD suite), `#66695`
(Feagin high-order ODE solvers, +5866/-35), `#66576`, `#67359`, `#66922`,
`#66886`, `#66880`, `#66850`, `#66842`, `#66834`, `#66742`, `#66767`,
`#66775`, `#66949`, `#66816`, `#68367`, `#60457`.

Everything else from the earlier list is in `Applied` or `Rejected` above.

The `source` rows (whole replacement files rather than diffs) are mostly
the 2025 specfunc/randist series from one contributor: `#66816`,
`#66834`, `#66842`, `#66850`, `#66862`, `#66877`, `#66880`, `#66886`,
`#66922`, `#66949`, `#67359`, `#67728`, `#66767`, `#66775`, `#66800`,
`#67621`, `#42472`.  These need the patch to be reconstructed by
hand; `#47345`, `#60371`, `#67621`, `#42472` and `#47646` came out of this
group and are now applied.


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
