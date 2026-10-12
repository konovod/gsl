noinst_LTLIBRARIES = libgslblas.la

pkginclude_HEADERS = gsl_blas.h gsl_blas_types.h

AM_CPPFLAGS = -I$(top_srcdir)

libgslblas_la_SOURCES = blas.c 

TESTS = $(check_PROGRAMS)
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslblas.la ../cblas/libgslcblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la


