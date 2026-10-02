.. index:: Airy functions
.. index:: Ai(x)
.. index:: Bi(x)

The Airy functions :math:`Ai(x)` and :math:`Bi(x)` are defined by the
integral representations,

.. only:: not texinfo

   .. math::

      Ai(x) & = {1\over\pi} \int_0^\infty \cos(t^3/3 + xt ) \,dt \\
      Bi(x) & = {1\over\pi} \int_0^\infty (e^{-t^3/3 + xt} + \sin(t^3/3 + xt)) \,dt

.. only:: texinfo

   | Ai(x) = (1/\pi) \int_0^\infty \cos((1/3) t^3 + xt) dt
   | Bi(x) = (1/\pi) \int_0^\infty (e^(-(1/3) t^3 + xt) + \sin((1/3) t^3 + xt)) dt

For further information see Abramowitz & Stegun, Section 10.4. The Airy
functions are defined in the header file :file:`gsl_sf_airy.h`.

Airy Functions
--------------

.. function:: double gsl_sf_airy_Ai (double x, gsl_mode_t mode)
              int gsl_sf_airy_Ai_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the Airy function :math:`Ai(x)` with an accuracy
   specified by :data:`mode`.

   For large positive :data:`x` the result underflows and
   ``GSL_EUNDRFLW`` is returned.  For :data:`x < 0` the function
   oscillates, and once the phase can no longer be determined from a
   double the result is set to zero with the modulus as its error and
   ``GSL_ELOSS`` is returned; see the note below.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_ELOSS

.. function:: double gsl_sf_airy_Bi (double x, gsl_mode_t mode)
              int gsl_sf_airy_Bi_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the Airy function :math:`Bi(x)` with an accuracy
   specified by :data:`mode`.

   For large positive :data:`x` the result overflows and
   ``GSL_EOVRFLW`` is returned; for large negative :data:`x` the result
   is indeterminate and ``GSL_ELOSS`` is returned, as for
   :func:`gsl_sf_airy_Ai`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_EOVRFLW, GSL_ELOSS

.. function:: double gsl_sf_airy_Ai_scaled (double x, gsl_mode_t mode)
              int gsl_sf_airy_Ai_scaled_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute a scaled version of the Airy function
   :math:`S_A(x) Ai(x)`.  For :math:`x > 0` the scaling factor :math:`S_A(x)` is
   :math:`\exp(+(2/3) x^{3/2})`, 
   and is 1 for :math:`x < 0`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_ELOSS

.. function:: double gsl_sf_airy_Bi_scaled (double x, gsl_mode_t mode)
              int gsl_sf_airy_Bi_scaled_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute a scaled version of the Airy function
   :math:`S_B(x) Bi(x)`.  For :math:`x > 0` the scaling factor :math:`S_B(x)` is
   :math:`exp(-(2/3) x^{3/2})`, and is 1 for :math:`x < 0`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_ELOSS

.. note::

   For :math:`x < 0` the accuracy of the Airy functions and their
   derivatives is limited by the phase
   :math:`\theta \approx (2/3)|x|^{3/2}`.  Computing it in double
   precision leaves an absolute error of order
   :math:`\varepsilon\,|\theta|`, so the relative error of the result
   grows roughly like :math:`\varepsilon\,|x|^{3/2}`.  When that error
   reaches a radian the value is no longer determined by the argument and
   the routines return zero with the modulus (or, for a derivative, the
   amplitude) as the error, together with ``GSL_ELOSS``.  For
   :math:`x > 0` the accuracy is a few ulp.


Derivatives of Airy Functions
-----------------------------

.. function:: double gsl_sf_airy_Ai_deriv (double x, gsl_mode_t mode)
              int gsl_sf_airy_Ai_deriv_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the Airy function derivative :math:`Ai'(x)` with
   an accuracy specified by :data:`mode`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_ELOSS

.. function:: double gsl_sf_airy_Bi_deriv (double x, gsl_mode_t mode)
              int gsl_sf_airy_Bi_deriv_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the Airy function derivative :math:`Bi'(x)` with
   an accuracy specified by :data:`mode`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_EOVRFLW, GSL_ELOSS

.. function:: double gsl_sf_airy_Ai_deriv_scaled (double x, gsl_mode_t mode)
              int gsl_sf_airy_Ai_deriv_scaled_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the scaled Airy function derivative 
   :math:`S_A(x) Ai'(x)`.  
   For :math:`x > 0` the scaling factor :math:`S_A(x)` is
   :math:`\exp(+(2/3) x^{3/2})`, and is 1 for :math:`x < 0`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_ELOSS

.. function:: double gsl_sf_airy_Bi_deriv_scaled (double x, gsl_mode_t mode)
              int gsl_sf_airy_Bi_deriv_scaled_e (double x, gsl_mode_t mode, gsl_sf_result * result)

   These routines compute the scaled Airy function derivative 
   :math:`S_B(x) Bi'(x)`.
   For :math:`x > 0` the scaling factor :math:`S_B(x)` is
   :math:`exp(-(2/3) x^{3/2})`, and is 1 for :math:`x < 0`.
.. Domain: -infinity < x < infinity
.. Exceptional Return Values: GSL_ELOSS

Zeros of Airy Functions
-----------------------

.. function:: double gsl_sf_airy_zero_Ai (unsigned int s)
              int gsl_sf_airy_zero_Ai_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th zero of the Airy
   function :math:`Ai(x)`.

.. function:: double gsl_sf_airy_zero_Bi (unsigned int s)
              int gsl_sf_airy_zero_Bi_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th zero of the Airy
   function :math:`Bi(x)`.

Zeros of Derivatives of Airy Functions
--------------------------------------

.. function:: double gsl_sf_airy_zero_Ai_deriv (unsigned int s)
              int gsl_sf_airy_zero_Ai_deriv_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th zero of the Airy
   function derivative :math:`Ai'(x)`.

.. function:: double gsl_sf_airy_zero_Bi_deriv (unsigned int s)
              int gsl_sf_airy_zero_Bi_deriv_e (unsigned int s, gsl_sf_result * result)

   These routines compute the location of the :data:`s`-th zero of the Airy
   function derivative :math:`Bi'(x)`.
