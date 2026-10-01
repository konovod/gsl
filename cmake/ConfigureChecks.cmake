# ConfigureChecks.cmake - feature detection for the CMake build of GSL.
#
# Sets the HAVE_* variables consumed by cmake/config.h.cmake.  This mirrors
# the checks performed by configure.ac as closely as is practical without a
# shell / autoconf.

include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckCSourceCompiles)
include(CheckCSourceRuns)

# ---------------------------------------------------------------------------
# Headers
# ---------------------------------------------------------------------------
check_include_file(ieeefp.h HAVE_IEEEFP_H)
check_include_file(complex.h HAVE_COMPLEX_H)
check_include_file(unistd.h HAVE_UNISTD_H)

# STDC_HEADERS: the standard C headers are present (autoconf's AC_HEADER_STDC).
# Every C99+ toolchain we target satisfies this, and results.c relies on it to
# choose <stdarg.h> over the obsolete <varargs.h>.
set(STDC_HEADERS 1)

# ---------------------------------------------------------------------------
# Declared functions
# ---------------------------------------------------------------------------

# On glibc systems the declarations of hypot/log1p/expm1/acosh/asinh/atanh
# (and others) are only visible with a GNU/POSIX feature test macro.  GSL's
# own sources set _GNU_SOURCE via utils/system.h, and autoconf performs its
# probes in GNU mode, so use the same setting here.
if(NOT MSVC)
  add_compile_definitions(_GNU_SOURCE)
  set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE")
  # The math declarations live in libm on many systems; the linking test in
  # check_symbol_exists needs it or every probe fails at link time.
  include(CheckLibraryExists)
  check_library_exists(m cos "" GSL_HAVE_LIBM)
  if(GSL_HAVE_LIBM)
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif()
endif()

if(MSVC)
  # <fenv.h> is C99 and not provided by MSVC, so neither function exists.
  set(HAVE_DECL_FEENABLEEXCEPT 0)
  set(HAVE_DECL_FESETTRAPENABLE 0)
else()
  check_symbol_exists(feenableexcept  "fenv.h" HAVE_DECL_FEENABLEEXCEPT)
  check_symbol_exists(fesettrapenable "fenv.h" HAVE_DECL_FESETTRAPENABLE)
endif()

check_symbol_exists(vprintf stdio.h HAVE_VPRINTF)

# ---------------------------------------------------------------------------
# Declared math functions
# ---------------------------------------------------------------------------
check_symbol_exists(hypot    math.h HAVE_DECL_HYPOT)
check_symbol_exists(log1p    math.h HAVE_DECL_LOG1P)
check_symbol_exists(expm1    math.h HAVE_DECL_EXPM1)
check_symbol_exists(acosh    math.h HAVE_DECL_ACOSH)
check_symbol_exists(asinh    math.h HAVE_DECL_ASINH)
check_symbol_exists(atanh    math.h HAVE_DECL_ATANH)
check_symbol_exists(ldexp    math.h HAVE_DECL_LDEXP)
check_symbol_exists(frexp    math.h HAVE_DECL_FREXP)
check_symbol_exists(isinf    math.h HAVE_DECL_ISINF)
check_symbol_exists(isfinite math.h HAVE_DECL_ISFINITE)
check_symbol_exists(isnan    math.h HAVE_DECL_ISNAN)

if(HAVE_IEEEFP_H)
  check_symbol_exists(finite "math.h;ieeefp.h" HAVE_DECL_FINITE)
else()
  check_symbol_exists(finite math.h HAVE_DECL_FINITE)
endif()

check_c_source_compiles("
#include <stdlib.h>
int main(void) { return EXIT_SUCCESS; }
" HAVE_EXIT_SUCCESS_AND_FAILURE)

# ---------------------------------------------------------------------------
# printf long double support (runtime check when not cross compiling)
# ---------------------------------------------------------------------------
if(CMAKE_CROSSCOMPILING)
  set(HAVE_PRINTF_LONGDOUBLE 1)
else()
  check_c_source_runs("
#include <stdio.h>
int main(void)
{
  long double x = 1.234;
  if (printf(\"%Lg\", x) < 0) return 1;
  return 0;
}
" HAVE_PRINTF_LONGDOUBLE)
endif()

# ---------------------------------------------------------------------------
# inline semantics
#
# MSVC has no source-compatible GNU/C99 inline, so the *\/inline.c fallback
# translation units are compiled instead and HAVE_INLINE stays undefined.
# ---------------------------------------------------------------------------
if(MSVC)
  set(HAVE_INLINE 0)
else()
  check_c_source_compiles("
inline int gsl_cmake_inline_probe(void) { return 0; }
int main(void) { return gsl_cmake_inline_probe(); }
" HAVE_INLINE)
  if(HAVE_INLINE)
    check_c_source_compiles("
inline int gsl_cmake_inline_probe(void) { return 0; }
" GSL_C99_INLINE_PROBE)
    if(GSL_C99_INLINE_PROBE)
      set(HAVE_C99_INLINE 1)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------
# Extended precision registers
# ---------------------------------------------------------------------------
if(MSVC)
  set(HAVE_EXTENDED_PRECISION_REGISTERS 0)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i[3-6]86|x86_64|AMD64|amd64")
  set(HAVE_EXTENDED_PRECISION_REGISTERS 1)
else()
  set(HAVE_EXTENDED_PRECISION_REGISTERS 0)
endif()

# ---------------------------------------------------------------------------
# IEEE arithmetic interface
#
# The selected interface chooses which ieee-utils/fp-*.c is compiled.  MSVC
# has no interface file; ieee-utils/env.c uses cmake/msvc-compat/fp-win.c instead.
# ---------------------------------------------------------------------------
set(GSL_IEEE_INTERFACE "")
if(MSVC)
  set(GSL_IEEE_INTERFACE "")
elseif(APPLE)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
    set(GSL_IEEE_INTERFACE "")
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|i386|i686")
    set(GSL_IEEE_INTERFACE darwin86)
  else()
    set(GSL_IEEE_INTERFACE darwin)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "i[3-6]86|x86_64|AMD64|amd64")
    set(GSL_IEEE_INTERFACE gnux86)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "sparc")
    set(GSL_IEEE_INTERFACE gnusparc)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc|powerpc")
    set(GSL_IEEE_INTERFACE gnuppc)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "m68k")
    set(GSL_IEEE_INTERFACE gnum68k)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
  set(GSL_IEEE_INTERFACE freebsd)
elseif(CMAKE_SYSTEM_NAME STREQUAL "NetBSD")
  set(GSL_IEEE_INTERFACE netbsd)
elseif(CMAKE_SYSTEM_NAME STREQUAL "OpenBSD")
  set(GSL_IEEE_INTERFACE openbsd)
elseif(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
  set(GSL_IEEE_INTERFACE solaris)
elseif(CMAKE_SYSTEM_NAME STREQUAL "AIX")
  set(GSL_IEEE_INTERFACE aix)
endif()

if(GSL_IEEE_INTERFACE)
  string(TOUPPER "${GSL_IEEE_INTERFACE}" _gsl_iface_upper)
  set(HAVE_${_gsl_iface_upper}_IEEE_INTERFACE 1)
  message(STATUS "GSL IEEE arithmetic interface: ${GSL_IEEE_INTERFACE}")
else()
  message(STATUS "GSL IEEE arithmetic interface: none")
endif()

# ---------------------------------------------------------------------------
# IEEE comparisons (used as a fallback in sys/infnan.c)
# ---------------------------------------------------------------------------
if(CMAKE_CROSSCOMPILING)
  set(HAVE_IEEE_COMPARISONS 1)
else()
  check_c_source_runs("
int main(void)
{
  double x = 0.0;
  double nan = x / x;   /* generate a NaN without a constant expression */
  return (nan != nan) ? 0 : 1;
}
" HAVE_IEEE_COMPARISONS)
endif()
