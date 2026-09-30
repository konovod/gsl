/* config.h.cmake - template for the CMake build of GSL.
 *
 * The fragments below mirror the parts of the autoconf-generated config.h
 * that the library sources depend on (the AH_BOTTOM content from
 * configure.ac plus the AC_DEFINE'd feature macros).
 */

#ifndef GSL_CMAKE_CONFIG_H
#define GSL_CMAKE_CONFIG_H

/* --- header checks ------------------------------------------------------- */
#cmakedefine HAVE_IEEEFP_H 1
#cmakedefine HAVE_COMPLEX_H 1
#cmakedefine HAVE_UNISTD_H 1
#cmakedefine STDC_HEADERS 1

/* --- declared functions --------------------------------------------------- */
#cmakedefine HAVE_VPRINTF 1
#cmakedefine HAVE_DECL_FEENABLEEXCEPT 1
#cmakedefine HAVE_DECL_FESETTRAPENABLE 1

/* --- declared math functions --------------------------------------------- */
#cmakedefine HAVE_DECL_HYPOT 1
#cmakedefine HAVE_DECL_LOG1P 1
#cmakedefine HAVE_DECL_EXPM1 1
#cmakedefine HAVE_DECL_ACOSH 1
#cmakedefine HAVE_DECL_ASINH 1
#cmakedefine HAVE_DECL_ATANH 1
#cmakedefine HAVE_DECL_LDEXP 1
#cmakedefine HAVE_DECL_FREXP 1
#cmakedefine HAVE_DECL_ISINF 1
#cmakedefine HAVE_DECL_ISFINITE 1
#cmakedefine HAVE_DECL_FINITE 1
#cmakedefine HAVE_DECL_ISNAN 1

/* --- compiler / runtime characteristics ---------------------------------- */
#cmakedefine HAVE_INLINE 1
#cmakedefine HAVE_C99_INLINE 1
#cmakedefine HAVE_PRINTF_LONGDOUBLE 1
#cmakedefine HAVE_EXTENDED_PRECISION_REGISTERS 1
#cmakedefine HAVE_IEEE_COMPARISONS 1
#cmakedefine HAVE_IEEE_DENORMALS 1
#cmakedefine HAVE_FPU_X86_SSE 1

/* --- selected IEEE arithmetic interface (at most one is defined) --------- */
#cmakedefine HAVE_GNUSPARC_IEEE_INTERFACE 1
#cmakedefine HAVE_GNUM68K_IEEE_INTERFACE 1
#cmakedefine HAVE_GNUPPC_IEEE_INTERFACE 1
#cmakedefine HAVE_GNUX86_IEEE_INTERFACE 1
#cmakedefine HAVE_SUNOS4_IEEE_INTERFACE 1
#cmakedefine HAVE_SOLARIS_IEEE_INTERFACE 1
#cmakedefine HAVE_HPUX11_IEEE_INTERFACE 1
#cmakedefine HAVE_HPUX_IEEE_INTERFACE 1
#cmakedefine HAVE_TRU64_IEEE_INTERFACE 1
#cmakedefine HAVE_IRIX_IEEE_INTERFACE 1
#cmakedefine HAVE_AIX_IEEE_INTERFACE 1
#cmakedefine HAVE_FREEBSD_IEEE_INTERFACE 1
#cmakedefine HAVE_OS2EMX_IEEE_INTERFACE 1
#cmakedefine HAVE_NETBSD_IEEE_INTERFACE 1
#cmakedefine HAVE_OPENBSD_IEEE_INTERFACE 1
#cmakedefine HAVE_DARWIN_IEEE_INTERFACE 1
#cmakedefine HAVE_DARWIN86_IEEE_INTERFACE 1

/* --- AH_BOTTOM content from configure.ac --------------------------------- */

/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#if !HAVE_EXIT_SUCCESS_AND_FAILURE
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

/* Define a rounding function which moves extended precision values out of
   registers and rounds them to double precision. */
#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

/* Substitute gsl functions for missing system functions. */
#if !HAVE_DECL_HYPOT
#define hypot gsl_hypot
#endif

#if !HAVE_DECL_LOG1P
#define log1p gsl_log1p
#endif

#if !HAVE_DECL_EXPM1
#define expm1 gsl_expm1
#endif

#if !HAVE_DECL_ACOSH
#define acosh gsl_acosh
#endif

#if !HAVE_DECL_ASINH
#define asinh gsl_asinh
#endif

#if !HAVE_DECL_ATANH
#define atanh gsl_atanh
#endif

#if !HAVE_DECL_LDEXP
#define ldexp gsl_ldexp
#endif

#if !HAVE_DECL_FREXP
#define frexp gsl_frexp
#endif

#if !HAVE_DECL_ISINF
#define isinf gsl_isinf
#endif

#if !HAVE_DECL_ISFINITE
#define isfinite gsl_finite
#endif

#if !HAVE_DECL_FINITE
#define finite gsl_finite
#endif

#if !HAVE_DECL_ISNAN
#define isnan gsl_isnan
#endif

#ifdef __GNUC__
#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);
#else
#define DISCARD_POINTER(p) /* ignoring discarded pointer */
#endif

#if defined(GSL_RANGE_CHECK_OFF) || !defined(GSL_RANGE_CHECK)
#define GSL_RANGE_CHECK 0  /* turn off range checking by default internally */
#endif

#define RETURN_IF_NULL(x) if (!x) { return ; }

#endif /* GSL_CMAKE_CONFIG_H */
