/* getopt.h - minimal compatibility shim for the MSVC build.
 *
 * utils/system.h includes <getopt.h> when HAVE_UNISTD_H is defined.
 * Nothing in libgsl actually uses getopt, so declaring the standard
 * entry point is enough to let system.h compile.
 */

#ifndef GSL_MSVC_GETOPT_H
#define GSL_MSVC_GETOPT_H

#ifdef __cplusplus
extern "C" {
#endif

extern char *optarg;
extern int optind;
extern int opterr;
extern int optopt;

int getopt (int argc, char *const *argv, const char *shortopts);

#ifdef __cplusplus
}
#endif

#endif /* GSL_MSVC_GETOPT_H */
