GSL @GSL_VERSION@ - Windows x64 binaries
Built from commit @GSL_SHORT_SHA@ with MSVC, static CRT.

Layout
  include/gsl/*.h    public headers
  lib/gsl.lib        import library for gsl.dll
  lib/gslcblas.lib   import library for gslcblas.dll
  bin/gsl.dll, bin/gslcblas.dll
  bin/gsl-randist.exe, bin/gsl-histogram.exe

Usage
  cl /I<prefix>/include /DWIN32 /DGSL_DLL app.c
     /link /LIBPATH:<prefix>/lib gsl.lib gslcblas.lib
  and put <prefix>/bin on PATH, or copy the DLLs next to your executable.

  Link against the .lib import libraries, not the DLLs.

  /DGSL_DLL makes GSL's data symbols (GSL_VAR) import from the DLL instead of
  being declared extern; it is required when linking against the import
  library.

Notes
  The DLLs link the static MSVC runtime, so the Visual C++ redistributable is
  not required.  The import libraries are ordinary COFF objects and can also be
  consumed from MinGW or clang-cl.