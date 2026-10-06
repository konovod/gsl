.. index:: Zeta functions

The Riemann zeta function is defined in Abramowitz & Stegun, Section
23.2.  The functions described in this section are declared in the
header file :file:`gsl_sf_zeta.h`.

Riemann Zeta Function
---------------------
.. index:: Riemann Zeta Function

The Riemann zeta function is defined by the infinite sum

.. math:: \zeta(s) = \sum_{k=1}^\infty k^{-s}

.. function:: double gsl_sf_zeta_int (int n)
              int gsl_sf_zeta_int_e (int n, gsl_sf_result * result)

   These routines compute the Riemann zeta function :math:`\zeta(n)` 
   for integer :data:`n`,
   :math:`n \ne 1`.
.. Domain: n integer, n != 1
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_zeta (double s)
              int gsl_sf_zeta_e (double s, gsl_sf_result * result)

   These routines compute the Riemann zeta function :math:`\zeta(s)`
   for arbitrary :data:`s`,
   :math:`s \ne 1`.
.. Domain: s != 1.0
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: int gsl_sf_complex_zeta_e (double x, double y, gsl_sf_result * result_re, gsl_sf_result * result_im)

   This function computes the Riemann zeta function :math:`\zeta(s)` for the
   complex argument :math:`s = x + iy`, storing the real and imaginary parts in
   :data:`result_re` and :data:`result_im`.  The trivial zeros at the negative
   even integers are returned exactly.  The value is computed with the
   Euler-Maclaurin summation formula for :math:`\Re(s) \ge 0` and with the
   reflection formula :math:`\zeta(s) = 2^s \pi^{s-1} \sin(\pi s/2)
   \Gamma(1-s) \zeta(1-s)` for :math:`\Re(s) < 0`.
.. Domain: s != 1.0
.. Exceptional Return Values: GSL_EDOM

Riemann Zeta Function Minus One
-------------------------------

For large positive argument, the Riemann zeta function approaches one.
In this region the fractional part is interesting, and therefore we
need a function to evaluate it explicitly.

.. function:: double gsl_sf_zetam1_int (int n)
              int gsl_sf_zetam1_int_e (int n, gsl_sf_result * result)

   These routines compute :math:`\zeta(n) - 1` for integer :data:`n`,
   :math:`n \ne 1`.
.. Domain: n integer, n != 1
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

.. function:: double gsl_sf_zetam1 (double s)
              int gsl_sf_zetam1_e (double s, gsl_sf_result * result)

   These routines compute :math:`\zeta(s) - 1` for arbitrary :data:`s`,
   :math:`s \ne 1`.
.. Domain: s != 1.0
.. Exceptional Return Values: GSL_EDOM, GSL_EOVRFLW

Hurwitz Zeta Function
---------------------
.. index:: Hurwitz Zeta Function

The Hurwitz zeta function is defined by

.. math:: \zeta(s,q) = \sum_0^\infty (k+q)^{-s}

.. function:: double gsl_sf_hzeta (double s, double q)
              int gsl_sf_hzeta_e (double s, double q, gsl_sf_result * result)

   These routines compute the Hurwitz zeta function :math:`\zeta(s,q)` for
   :math:`s > 1`, :math:`q > 0`.
.. Domain: s > 1.0, q > 0.0
.. Exceptional Return Values: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: int gsl_sf_complex_hzeta_e (double x, double y, double q, gsl_sf_result * result_re, gsl_sf_result * result_im)

   This function computes the Hurwitz zeta function :math:`\zeta(s,q)` for the
   complex argument :math:`s = x + iy` and real :math:`q > 0`, storing the real
   and imaginary parts in :data:`result_re` and :data:`result_im`.  It uses the
   Euler-Maclaurin summation formula and is defined for :math:`\Re(s) > 1`.
.. Domain: Re(s) > 1.0, q > 0.0
.. Exceptional Return Values: GSL_EDOM

Eta Function
------------
.. index:: Eta Function

The eta function is defined by

.. math:: \eta(s) = (1-2^{1-s}) \zeta(s)

.. function:: double gsl_sf_eta_int (int n)
              int gsl_sf_eta_int_e (int n, gsl_sf_result * result)

   These routines compute the eta function :math:`\eta(n)` for integer :data:`n`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: double gsl_sf_eta (double s)
              int gsl_sf_eta_e (double s, gsl_sf_result * result)

   These routines compute the eta function :math:`\eta(s)` for arbitrary :data:`s`.
.. Exceptional Return Values: GSL_EUNDRFLW, GSL_EOVRFLW

.. function:: int gsl_sf_complex_eta_e (double x, double y, gsl_sf_result * result_re, gsl_sf_result * result_im)

   This function computes the eta function :math:`\eta(s)` for the complex
   argument :math:`s = x + iy` through :math:`\eta(s) = (1 - 2^{1-s})
   \zeta(s)`, storing the real and imaginary parts in :data:`result_re` and
   :data:`result_im`.  The value :math:`\eta(1) = \log 2` is handled exactly.
.. Domain: s != 1.0
.. Exceptional Return Values: GSL_EDOM
