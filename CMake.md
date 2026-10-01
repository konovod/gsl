Building GSL with CMake
=======================

This directory tree can also be built with CMake (>= 3.22), independently
of the autotools and Visual Studio project files, which remain available.

Supported platforms
-------------------

  * Windows with MSVC (tested with Visual Studio 2022 Build Tools, x64)
  * Linux (tested with GCC 13 and Clang, Ubuntu)
  * macOS (tested with Clang)

Quick start (Windows / MSVC)
----------------------------

From a "x64 Native Tools Command Prompt", or by first running vcvars64.bat:

    cmake -S . -B build-cmake -G Ninja -DCMAKE_BUILD_TYPE=Release
    cmake --build build-cmake
    ctest --test-dir build-cmake --output-on-failure     (tests, see below)

This produces gsl.dll, gslcblas.dll, their import libraries, and the
gsl-randist / gsl-histogram programs.  `-G "Visual Studio 17 2022" -A x64`
works as well; then use `cmake --build build-cmake --config Release`.

For a static build add `-DBUILD_SHARED_LIBS=OFF`.

Options
-------

  BUILD_SHARED_LIBS      ON      build shared libraries (DLL) rather than static
  GSL_BUILD_TOOLS        ON      build gsl-randist and gsl-histogram
  GSL_BUILD_TESTS        OFF     build and register the per-module test programs
  GSL_ENABLE_RANGE_CHECK OFF     enable run-time range checking inside the library

Installing and consuming
------------------------

    cmake --install build-cmake --prefix C:/some/prefix

installs headers under include/gsl, the libraries, the command line tools,
gsl.pc, gsl-config and a CMake package.

CMake consumers should request the config package explicitly, so that
CMake's own legacy FindGSL.cmake module is not used instead:

    find_package(GSL CONFIG REQUIRED)
    target_link_libraries(app PRIVATE GSL::gsl GSL::gslcblas)

The imported targets propagate the `WIN32;GSL_DLL` preprocessor definitions
needed for correct import of data symbols (GSL_VAR) when linking against
the DLL.

Tests
-----

`-DGSL_BUILD_TESTS=ON` builds one executable per automake check_PROGRAMS
entry and registers them with CTest.  Each test runs in its module source
directory because several tests read and write data files there, so run
them through CTest rather than invoking the binaries by hand.

Continuous integration
----------------------

.github/workflows/cmake.yml builds and tests the project on Linux (GCC and
Clang, shared and static), macOS (Clang, shared and static) and Windows
(MSVC, shared and static), and separately checks `cmake --install` followed
by a downstream `find_package(GSL CONFIG)` consumer on each platform.

Notes on the implementation
---------------------------
  * cmake/regen_sources.py regenerates cmake/gsl_sources.cmake from the
    automake *_la_SOURCES variables, and cmake/regen_tests.py regenerates
    cmake/gsl_tests.cmake from check_PROGRAMS.  Run them after changing a
    Makefile.am source list.

  * On MSVC the DLL exports are produced with WINDOWS_EXPORT_ALL_SYMBOLS;
    the gslhdrs / gsldefs header rewriting used by the Visual Studio
    projects is not needed.  ieee-utils uses cmake/msvc-compat/fp-win.c on
    Windows (there is no upstream Windows IEEE interface file).

  * cmake/msvc-compat provides minimal <unistd.h> and <getopt.h> shims for
    Windows; some test sources include them unconditionally.
