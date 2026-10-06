# AGENTS.md — working in this fork of GSL

Guidance for coding agents. Read this before making changes. It assumes the
rest of the repository is unfamiliar to you.

## What this repository is

This is a fork of **GNU GSL** (the GNU Scientific Library). It exists to
keep the library usable on current toolchains and to fold in the worthwhile
patches that have piled up in the Savannah bug tracker. Submitting work
upstream is desirable but **not** the immediate priority, so the fork is
maintained in a way that lets upstream-suitable changes be separated out
later.

Three kinds of work live side by side here:

1. **Patches harvested from Savannah** — bug fixes, documentation fixes and
   test improvements. Most are intended to go upstream eventually.
2. **Fork-only infrastructure** — the CMake build and GitHub Actions CI.
   Upstream maintains autotools; these are expected to stay in the fork.
3. **Rejected candidates** — patches that were evaluated and deliberately
   not applied, kept only as a record.

The classification of every change is recorded in `FORKNEWS`.

## Golden rules

1. **Every change is classified and recorded.** After implementing a change,
   add it to `FORKNEWS` with a disposition tag (`[upstream]`, `[fork]` or
   `[rejected]`). A change that is not in `FORKNEWS` does not exist as far as
   this fork is concerned.
2. **Keep upstream-bound and fork-only work in separate commits.** A
   `[fork]` change and an `[upstream]` change must never share a commit, even
   if they touch the same file. This is the whole point of the split, so it
   can be cherry-picked to `savannah/master` later.
3. **Update `FORKNEWS` in its own commit.** The convention is a commit
   titled `FORKNEWS: record <what>` immediately after the substantive
   commit(s). Do not amend the change commit.
4. **Never add a fork-only change to upstream's `NEWS`.** `NEWS` is
   upstream's release log; only touch it when a change is actually merged
   upstream.
5. **Never push to `origin` or `savannah`.** Only the `fork` remote is a push
   target, and only when explicitly asked.
6. **Verify before you believe.** `git apply --check` proves a patch still
   applies; it says nothing about correctness. Reproduce the behaviour with
   the built library and check against an independent reference.

## Remotes and branches

| remote | repository | role |
|---|---|---|
| `fork` | `https://github.com/konovod/gsl.git` | this fork; the push target |
| `origin` | `https://github.com/BrianGladman/gsl` | fork origin (CMake work) |
| `savannah` | `https://git.savannah.gnu.org/git/gsl.git` | upstream GNU GSL; read-only reference |

`master` tracks `fork/master` and is a direct descendant of
`savannah/master` with the fork's commits on top. Compare against
`savannah/master` when you need to know what is fork-specific.

## Fork bookkeeping: `FORKNEWS`

`FORKNEWS` (no extension) is the fork's change log, deliberately separate
from `NEWS`. It has two parts:

* **Summary** — one line per change:
  `date | tag | description | commit(s)`.
* **Per-change entries** — a short essay: what was wrong, what the fix does,
  which files/commits, how it was verified, and the Savannah bug id when
  there is one.

Tags and their meaning:

* `[upstream]` — intended for or suitable for upstream GSL; written in a
  style upstream would accept.
* `[fork]` — specific to this fork (CMake, CI) and not intended upstream.
* `[rejected]` — evaluated and deliberately not applied; the entry records
  the reasoning so a future reader does not re-investigate.

A `[fork]` change may be promoted to `[upstream]` later. When something is
actually merged upstream, remove its entry from `FORKNEWS` and add it to
`NEWS` if appropriate.

Commit subjects follow a `area: imperative summary` form, e.g.
`specfunc: fix the Hermite function derivative at orders 0 and 1`,
`linalg: initialize the LU pivot array before the factorization`,
`doc: document gsl_sf_complex_psi_e`. Bodies explain the defect, the fix and
the verification, wrapped at roughly 72–80 columns, and reference
`Savannah bug #NNNNN` (never GitHub `Fixes #NNN` syntax). Commits here carry
no `Signed-off-by` trailer.

## Savannah bug review

The tracker is <https://savannah.gnu.org/bugs/?group=gsl>.
Texts for all bugs are available offline at `/temp/savannah-store/dossiers/NNN.txt` 

* `scripts/savannah_bugs.py` walks the open bugs, downloads patch-like
  attachments, and runs `git apply --check` against the tree. It writes an
  inventory and the patches under a temp directory (default
  `%TEMP%/savannah-gsl`), never into the tree. A full sweep is a few minutes;
  results are cached in `cache.json`.
* `SAVANNAH_REVIEW.md` holds the working notes and verdicts for each bug.
  It is deliberately **not** part of the fork's record of changes and is
  kept out of `FORKNEWS`; it records both applied and rejected items, plus
  method notes.
* `/temp/` (git-ignored) is the scratch area for patch and reference files
  being evaluated.
* `SAVANNAH_TRIAGE.md` consists of two parts - big table with short description and status of each bug and backlog with items remaining to process.

**Eligibility rule.** A candidate is in scope if it is a **bug fix,
documentation correction, or test-quality improvement**. Feature requests,
new API, new algorithms and performance-only changes was out of scope and
deferred, however cleanly they apply. We are widening our eligibility rule now, starting from changes with minimal impact on backward compatibility.


**Review method** (this is the hard-won part — follow it):

1. Reproduce the reported behaviour with the built library before changing
   anything. Confirm the diagnosis rather than trusting the report: several
   Savannah patches apply cleanly and are wrong or actively harmful.
2. Reconstruct "whole replacement file" attachments by hand; there is no
   diff to apply.
3. Implement the fix in the fork's style, not necessarily verbatim.
4. Add a regression test that **fails with the fix reverted** (the negative
   control). Record the observed failure in `FORKNEWS`.
5. Check against an independent reference (a different algorithm, higher
   precision, or a published value), not against GSL itself.
6. Record the outcome: applied entries in the `FORKNEWS` summary and body;
   rejected ones as `[rejected]` with the reason, and in
   `SAVANNAH_REVIEW.md`. Update status and backlog in `SAVANNAH_TRIAGE.md`

## Build and test

### CMake (the fork's primary build)

Requires CMake >= 3.22. Out-of-source only; `build-cmake/` is git-ignored.

Linux / macOS:

```sh
cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release -DGSL_BUILD_TESTS=ON
cmake --build build-cmake --parallel
ctest --test-dir build-cmake --output-on-failure
```

Windows / MSVC — from an x64 Native Tools prompt, or after `vcvars64.bat`:

```sh
cmake -S . -B build-cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DGSL_BUILD_TESTS=ON
cmake --build build-cmake --parallel
ctest --test-dir build-cmake --output-on-failure
```

`cmake --build` and `ctest` on a Visual Studio generator build need
`--config Release` and `-C Release` respectively. An existing MSVC-configured
tree (e.g. `build-cmake/`) is rebuilt with:

```
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && cd /d <repo>\build-cmake && ninja'
```

Options (from `CMake.md`):

| option | default | meaning |
|---|---|---|
| `BUILD_SHARED_LIBS` | `ON` | shared library (DLL) vs static |
| `GSL_BUILD_TOOLS` | `ON` | build `gsl-randist` and `gsl-histogram` |
| `GSL_BUILD_TESTS` | `OFF` | build and register the per-module test programs |
| `GSL_ENABLE_RANGE_CHECK` | `OFF` | run-time range checking inside the library |

Run tests through CTest, not the executables directly: several tests read and
write data files in their module source directory, which CTest sets up.

Baseline expectation: **all `ctest` tests pass** (56/56 at the time of
writing). A change is not done until the full suite is green.

### Autotools

The upstream autotools build is still present and must not be broken. Avoid
changing `configure.ac`/`Makefile.am` unless the change is genuinely
upstream-bound, and keep it working when you do.

### Adding or removing a source file

A new library source must be added to **both** build systems:

1. the module's `Makefile.am` (`*_la_SOURCES` or `check_PROGRAMS`), and
2. `cmake/gsl_sources.cmake` / `cmake/gsl_tests.cmake`.

Do not edit the generated CMake lists by hand — regenerate them from the
automake variables:

```sh
python cmake/regen_sources.py
python cmake/regen_tests.py
```

The generated files say so in their header.

## Coding conventions

GSL is C89-flavoured portable C with GNU style. Follow the surrounding code;
do not introduce a new style.

* 2-space indentation, braces on their own line indented one level for
  functions and blocks; GNU brace placement.
* Prefer `/* ... */` comments; the codebase is C89-compatible.
* Every source starts with `#include <config.h>` before the system and GSL
  headers.
* Public symbols are prefixed `gsl_` / `GSL_`; internal helpers are
  `static`. Keep new internal names static unless exporting is required.
* Report errors with the `GSL_ERROR` / `GSL_ERROR_VAL` / `DOMAIN_ERROR`
  family and return the documented `GSL_*` code. Do not `printf` and do not
  call `abort`.
* **Do not add `isnan` guards.** In `specfunc` and elsewhere, a `NaN`
  argument is intentionally propagated, not treated as a domain error. A
  non-finite argument that has no finite limit (e.g. `sin_pi(inf)`) is a
  domain error (`GSL_EDOM`); one that does have a limit is not.
* Keep error estimates (`result->err`) meaningful. A reported error of
  exactly zero for an inexact computation is a bug — the test harness
  rejects expected values that fall outside the reported error bar.

## Tests

* Tests live next to the module (`<module>/test.c` and helpers) and use the
  `TEST_SF` / `TEST_*` harness with tolerances `TEST_TOL0` (tightest) through
  `TEST_TOL6`. Pick a tolerance that reflects the conditioning of the
  computation and the accuracy of the reference — not the tightest that
  happens to pass on your machine.
* Add a vector that specifically exercises the bug, and verify it fails with
  the fix reverted. Record the failure count in `FORKNEWS`.
* Expected values must lie inside the library's own reported error bar, or
  `TEST_SF` reports `TEST_SF_INCONS`. At exact zeros there is no relative
  scale; use an absolute-aware tolerance rather than `TEST_TOL0`.

## Documentation

There are two manual trees and they mirror each other:

* `doc/*.rst` — the Sphinx manual (current).
* `doc_texinfo/*.texi` — the Texinfo manual.

When a change alters documented behaviour, update **both** where the
corresponding entry exists, and update the examples under `doc/examples/`
and the stale mirror `doc_texinfo/examples/` together. Deprecations use the
`.. deprecated::` directive in Sphinx and equivalent prose in Texinfo, and
must be derived from the actual `GSL_DISABLE_DEPRECATED` blocks in the
headers — do not mark a function deprecated just because a comment mentions
it.

## Platforms and numerical pitfalls

This fork is built and tested on Windows (MSVC), Linux and macOS. Changes
must pass on all of them.

* **`long` width differs.** Windows/MSVC is LLP64 (`long` = 32 bit); Linux
  and macOS are LP64. Code that relies on `long` overflow or on
  `unsigned long` holding a 64-bit value is wrong on one family. Use
  `int32_t`/`uint32_t` or the `GSL_*` types.
* **Extended precision is not uniform.** x86-64 keeps intermediates in 80-bit
  x87 registers; MSVC and arm64 macOS round to `double`. A test that passes
  only because of x87 excess precision will fail on macOS; widen the
  tolerance or restructure the computation.
* **FMA / `-ffp-contract=fast`** is a default in many compiler setups and
  changes rounding. Tests must not depend on a single contraction mode.
  `-ffast-math` is unsupported by GSL and is never an acceptable fix.
* **Non-ASCII build paths.** The Windows temp path here contains non-ASCII
  characters, which can break linking test programs against the DLL import
  library. The reliable workaround (archive the `.obj` files with `ar`) is
  written up in the "Method notes" section of `SAVANNAH_REVIEW.md`.
* On MSVC the DLL exports are produced with `WINDOWS_EXPORT_ALL_SYMBOLS`;
  `cmake/msvc-compat/` supplies minimal `unistd.h`/`getopt.h` shims and the
  Windows IEEE implementation.

## Repository layout

```
<module>/          one directory per GSL module (min, specfunc, linalg, ...):
                   sources, public gsl_*.h, test.c, Makefile.am
doc/               Sphinx manual sources + examples
doc_texinfo/       Texinfo manual sources + examples (mirror)
cmake/             CMake modules, generated source/test lists, regen scripts
.github/workflows/ CMake and Windows DLL CI
scripts/           upstream helper scripts + savannah_bugs.py
FORKNEWS           fork change log (classified)
CMake.md           CMake build documentation
SAVANNAH_REVIEW.md bug-review working notes (scratch)
HACKING            upstream developer notes (releases, checks, portability)
```

## Do not

* Do not mix `[upstream]` and `[fork]` changes in one commit.
* Do not touch `NEWS` for fork-only changes.
* Do not edit `cmake/gsl_sources.cmake` or `cmake/gsl_tests.cmake` by hand.
* Do not take feature, new-API or performance patches from Savannah.
* Do not trust a patch because it applies cleanly; reproduce and test it.
* Do not push to `origin` or `savannah`.
* Do not delete or overwrite unfamiliar files or uncommitted work — check
  first.
