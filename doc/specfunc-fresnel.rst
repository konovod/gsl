.. index:: Fresnel integrals

The Fresnel integrals are defined by

.. math::

   C(x) = \int_0^x \cos(\pi t^2/2) \, dt, \qquad
   S(x) = \int_0^x \sin(\pi t^2/2) \, dt

Both functions are odd, :math:`C(-x) = -C(x)` and
:math:`S(-x) = -S(x)`.  The routines are declared in the header file
:file:`gsl_sf_fresnel.h`.

.. function:: double gsl_sf_fresnel_c (double x)
              int gsl_sf_fresnel_c_e (double x, gsl_sf_result * result)

   These routines compute the Fresnel cosine integral :math:`C(x)` for
   :data:`x`.

.. function:: double gsl_sf_fresnel_s (double x)
              int gsl_sf_fresnel_s_e (double x, gsl_sf_result * result)

   These routines compute the Fresnel sine integral :math:`S(x)` for
   :data:`x`.

The values are computed in full double precision for moderate arguments.
For large :math:`|x|` the phase :math:`\pi x^2/2` can no longer be
represented exactly and the reported error grows with :math:`|x|`.  An
infinite argument returns the error :macro:`GSL_EDOM`, since the
integrals have no limit there.
