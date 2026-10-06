.. index::
   single: Jacobi elliptic functions
   single: elliptic functions (Jacobi)

The Jacobian Elliptic functions are defined in Abramowitz & Stegun,
Chapter 16.  The functions are declared in the header file
:file:`gsl_sf_elljac.h`.

.. function:: int gsl_sf_elljac_e (double u, double m, double * sn, double * cn, double * dn)

   This function computes the Jacobian elliptic functions :math:`sn(u|m)`,
   :math:`cn(u|m)`, :math:`dn(u|m)` by descending Landen
   transformations.
.. Exceptional Return Values: GSL_EDOM

.. function:: int gsl_sf_elljac_arcsn_e (double sn, double m, gsl_sf_result * result)
              int gsl_sf_elljac_arccn_e (double cn, double m, gsl_sf_result * result)
              int gsl_sf_elljac_arcdn_e (double dn, double m, gsl_sf_result * result)

   These functions compute the principal value :math:`u` of the inverse
   Jacobian elliptic functions with the parameter :data:`m`
   (:math:`|m| \le 1`), so that :math:`sn(u|m) = sn`,
   :math:`cn(u|m) = cn` and :math:`dn(u|m) = dn` respectively.  The
   principal value satisfies :math:`0 \le u \le 2K(m)`.  The parameter
   :data:`m` is used rather than the modulus :math:`k = \sqrt{m}`, matching
   :func:`gsl_sf_elljac_e`; negative values of :data:`m` are supported.
   The inverse of :math:`dn` is not defined for :math:`m = 0`, where
   :math:`dn(u|0) = 1` for every :math:`u`.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_elljac_arcsn (double sn, double m)
              double gsl_sf_elljac_arccn (double cn, double m)
              double gsl_sf_elljac_arcdn (double dn, double m)

   These functions return the same values as the corresponding
   :code:`_e` functions, using the default error handler.
