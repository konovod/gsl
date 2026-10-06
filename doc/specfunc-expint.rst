.. index::
   single: exponential integrals
   single: integrals, exponential

Information on the exponential integrals can be found in Abramowitz &
Stegun, Chapter 5.  These functions are declared in the header file
:file:`gsl_sf_expint.h`.

Exponential Integral
--------------------
.. index:: E1(x), E2(x), Ei(x)

.. function:: double gsl_sf_expint_E1 (double x)
              int gsl_sf_expint_E1_e (double x, gsl_sf_result * result)

   These routines compute the exponential integral :math:`E_1(x)`,

   .. math:: E_1(x) := \Re \int_1^\infty dt \exp(-xt)/t.

   The function diverges to :math:`+\infty` as :math:`x \to 0`
   (from either side), and :math:`E_1(0)` is returned as
   :data:`GSL_POSINF` with zero error.
.. Domain: x != 0.0
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_expint_E2 (double x)
              int gsl_sf_expint_E2_e (double x, gsl_sf_result * result)

   These routines compute the second-order exponential integral :math:`E_2(x)`,

   .. math:: E_2(x) := \Re \int_1^\infty dt \exp(-xt)/t^2

.. Domain: x != 0.0
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_expint_En (int n, double x)
              int gsl_sf_expint_En_e (int n, double x, gsl_sf_result * result)

   These routines compute the exponential integral :math:`E_n(x)` of order :data:`n`, 

   .. math:: E_n(x) := \Re \int_1^\infty dt \exp(-xt)/t^n.

   For :math:`n > 0` the function is finite at :math:`x = 0`.  For
   :math:`n = 0`, :math:`E_0(0)` diverges to :math:`+\infty` and is
   returned as :data:`GSL_POSINF` with zero error.
.. Domain: x != 0.0
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

Ei(x)
-----

.. function:: double gsl_sf_expint_Ei (double x)
              int gsl_sf_expint_Ei_e (double x, gsl_sf_result * result)

   These routines compute the exponential integral :math:`Ei(x)`,

   .. only:: not texinfo

      .. math:: \hbox{Ei}(x) = - PV \left( \int_{-x}^\infty dt \exp(-t)/t \right)

   .. only:: texinfo

      ::

         Ei(x) = - PV(\int_{-x}^\infty dt \exp(-t)/t)

   where :math:`PV` denotes the principal value of the integral.
   :math:`Ei(x)` diverges to :math:`-\infty` as :math:`x \to 0`, and
   :math:`Ei(0)` is returned as :data:`GSL_NEGINF` with zero error.
.. Domain: x != 0.0
.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

Complex Exponential Integrals
-----------------------------
.. index::
   single: exponential integrals, complex
   single: E1(z)
   single: En(z)
   single: Ei(z)

.. function:: int gsl_sf_complex_expint_E1_e (double zr, double zi, gsl_sf_result * re, gsl_sf_result * im)
              int gsl_sf_complex_expint_En_e (int n, double zr, double zi, gsl_sf_result * re, gsl_sf_result * im)

   These routines compute the exponential integral :math:`E_n(z)` of order
   :data:`n` for the complex argument :math:`z = zr + i zi`,

   .. math:: E_n(z) := \int_1^\infty dt \exp(-z t)/t^n

   continued analytically to the cut plane :math:`-\pi < \arg(z) \le \pi`,
   with the branch cut along the negative real axis approached from above.
   The order must satisfy :math:`n \ge 1`.
   :math:`z = 0` with :math:`n = 1` is a branch point and is a domain
   error; for :math:`n > 1` the value at :math:`z = 0` is
   :math:`1/(n-1)`.

   The real and imaginary parts are returned in :data:`re` and :data:`im`.
   A NaN component of :data:`z` propagates, while :math:`z \to +\infty`
   along the real axis has the limit zero.

.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW, GSL_ELOSS, GSL_EMAXITER, GSL_EFAILED

.. function:: int gsl_sf_complex_expint_E1_scaled_e (double zr, double zi, gsl_sf_result * re, gsl_sf_result * im)
              int gsl_sf_complex_expint_En_scaled_e (int n, double zr, double zi, gsl_sf_result * re, gsl_sf_result * im)

   These routines compute the scaled exponential integral
   :math:`\exp(z) E_n(z)`, which removes the exponential behaviour of
   :math:`E_n(z)` in both half-planes.

.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW, GSL_ELOSS, GSL_EMAXITER, GSL_EFAILED

.. function:: int gsl_sf_complex_expint_Ei_e (double zr, double zi, gsl_sf_result * re, gsl_sf_result * im)

   This routine computes the exponential integral :math:`Ei(z)` for the
   complex argument :math:`z = zr + i zi`,

   .. math:: \hbox{Ei}(z) = - PV \left( \int_{-z}^\infty dt \exp(-t)/t \right)

   with the same branch cut along the negative real axis approached from
   above, so that :math:`Ei(x)` is real for real :math:`x > 0`.  The
   origin is a branch point and is a domain error.

.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW, GSL_ELOSS, GSL_EMAXITER, GSL_EFAILED

The complex exponential integrals are computed with the algorithm of
D. E. Amos, "Computation of Exponential Integrals of a Complex Argument",
*ACM Transactions on Mathematical Software* 16(2), 169--177 (1990),
Algorithm 683.

Hyperbolic Integrals
--------------------
.. index::
   single: hyperbolic integrals
   single: Shi(x)
   single: Chi(x)

.. function:: double gsl_sf_Shi (double x)
              int gsl_sf_Shi_e (double x, gsl_sf_result * result)

   These routines compute the integral
   
   .. only:: not texinfo
      
      .. math:: \hbox{Shi}(x) = \int_0^x dt \sinh(t)/t

   .. only:: texinfo

      ::

         Shi(x) = \int_0^x dt \sinh(t)/t

.. Exceptional Return Values: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_Chi (double x)
              int gsl_sf_Chi_e (double x, gsl_sf_result * result)

   These routines compute the integral
   
   .. only:: not texinfo
      
      .. math:: \hbox{Chi}(x) := \Re \left[ \gamma_E + \log(x) + \int_0^x dt (\cosh(t)-1)/t \right]

   .. only:: texinfo

      ::

         Chi(x) := \Re[ \gamma_E + \log(x) + \int_0^x dt (\cosh(t)-1)/t ]

   where :math:`\gamma_E` is the Euler constant (available as the macro :macro:`M_EULER`).
.. Domain: x != 0.0
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW

Ei_3(x)
-------

.. function:: double gsl_sf_expint_3 (double x)
              int gsl_sf_expint_3_e (double x, gsl_sf_result * result)

   These routines compute the third-order exponential integral
   
   .. only:: not texinfo
      
      .. math:: {\rm Ei}_3(x) = \int_0^x dt \exp(-t^3)

   .. only:: texinfo

      ::

         Ei_3(x) = \int_0^x dt \exp(-t^3)`
         
   for :math:`x \ge 0`.

.. Exceptional Return Values: GSL_EDOM

Trigonometric Integrals
-----------------------
.. index::
   single: trigonometric integrals
   single: Si(x)
   single: Ci(x)

.. function:: double gsl_sf_Si (const double x)
              int gsl_sf_Si_e (double x, gsl_sf_result * result)

   These routines compute the Sine integral
   
   .. only:: not texinfo
      
      .. math:: \hbox{Si}(x) = \int_0^x dt \sin(t)/t

   .. only:: texinfo

      ::

         Si(x) = \int_0^x dt \sin(t)/t

.. Exceptional Return Values: none
 
.. function:: double gsl_sf_Ci (const double x)
              int gsl_sf_Ci_e (double x, gsl_sf_result * result)

   These routines compute the Cosine integral
   
   .. only:: not texinfo
      
      .. math:: \hbox{Ci}(x) = -\int_x^\infty dt \cos(t)/t

   .. only:: texinfo

      ::

         Ci(x) = -\int_x^\infty dt \cos(t)/t}
         
   for :math:`x > 0`.
   The limit as :math:`x \to +\infty` is returned for ``+inf``:
   :math:`\hbox{Si}(+\infty) = \pi/2` and
   :math:`\hbox{Ci}(+\infty) = 0`.  ``Ci`` has no limit at infinity
   for an argument other than :math:`+\infty`, so ``-inf`` is outside
   the domain; ``Si`` has the limit :math:`-\pi/2` at ``-inf`` and
   returns it.
.. Domain: x > 0.0
.. Exceptional Return Values: GSL_EDOM

Arctangent Integral
-------------------
.. index::
   single: arctangent integral

.. function:: double gsl_sf_atanint (double x)
              int gsl_sf_atanint_e (double x, gsl_sf_result * result)

   These routines compute the Arctangent integral, which is defined as
   
   .. only:: not texinfo
      
      .. math:: \hbox{AtanInt}(x) = \int_0^x dt \arctan(t)/t

   .. only:: texinfo

      ::

         AtanInt(x) = \int_0^x dt \arctan(t)/t

.. Domain: 
.. Exceptional Return Values: 
