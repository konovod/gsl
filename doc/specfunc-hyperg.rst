.. index::
   single: hypergeometric functions
   single: confluent hypergeometric functions

Hypergeometric functions are described in Abramowitz & Stegun, Chapters
13 and 15.  These functions are declared in the header file
:file:`gsl_sf_hyperg.h`.

.. function:: double gsl_sf_hyperg_0F1 (double c, double x)
              int gsl_sf_hyperg_0F1_e (double c, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_0F_1(c,x)

   .. only:: texinfo

      .. math:: 0F1(c,x)

.. It is related to Bessel functions
.. 0F1[c,x] =
..   Gamma[c]    x^(1/2(1-c)) I_(c-1)(2 Sqrt[x])
..   Gamma[c] (-x)^(1/2(1-c)) J_(c-1)(2 Sqrt[-x])
.. exceptions: GSL_EOVRFLW, GSL_EUNDRFLW

.. function:: double gsl_sf_hyperg_1F1_int (int m, int n, double x)
              int gsl_sf_hyperg_1F1_int_e (int m, int n, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(m,n,x) = M(m,n,x)

   .. only:: texinfo

      .. math:: 1F1(m,n,x) = M(m,n,x)

   for integer parameters :data:`m`, :data:`n`.

   The function is undefined when :data:`n` is zero or a negative
   integer, and these routines return a domain error in that case.  The
   exception is when the series terminates before the denominator
   vanishes, i.e. when :data:`m` is a nonpositive integer with
   :math:`m \ge n`; the value is then finite, and in particular
   :math:`M(m,n,0) = 1`.
.. exceptions: 

.. function:: double gsl_sf_hyperg_1F1 (double a, double b, double x)
              int gsl_sf_hyperg_1F1_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function

   .. only:: not texinfo

      .. math:: {}_1F_1(a,b,x) = M(a,b,x)

   .. only:: texinfo

      .. math:: 1F1(a,b,x) = M(a,b,x)

   for general parameters :data:`a`, :data:`b`.

   As with the integer-parameter form, the function is undefined when
   :data:`b` is zero or a negative integer, and these routines return a
   domain error in that case.  When the series terminates before the
   denominator vanishes, i.e. when :data:`a` is a nonpositive integer
   with :math:`a \ge b`, the value is finite, and in particular
   :math:`M(a,b,0) = 1`.

    A :code:`NaN` argument is propagated, as it is throughout this
    chapter; it is not treated as a domain error.
.. exceptions:

.. function:: double gsl_sf_hyperg_U_int (int m, int n, double x)
              int gsl_sf_hyperg_U_int_e (int m, int n, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n`.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_int_e10_e (int m, int n, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(m,n,x)` for integer parameters :data:`m`, :data:`n` using the
   :type:`gsl_sf_result_e10` type to return a result with extended range.

.. function:: double gsl_sf_hyperg_U (double a, double b, double x)
              int gsl_sf_hyperg_U_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the confluent hypergeometric function :math:`U(a,b,x)`.

   For :math:`x < 0` the principal value of :math:`U(a,b,x)` is complex,
   with a branch cut on the negative real axis.  These real-valued
   routines return the real part of the principal value, the real
   continuous solution of Kummer's equation.
.. exceptions:

.. function:: int gsl_sf_hyperg_U_e10_e (double a, double b, double x, gsl_sf_result_e10 * result)

   This routine computes the confluent hypergeometric function
   :math:`U(a,b,x)` using the :type:`gsl_sf_result_e10` type to return a
   result with extended range. 
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1 (double a, double b, double c, double x)
              int gsl_sf_hyperg_2F1_e (double a, double b, double c, double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) = F(a,b,c,x)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) = F(a,b,c,x)
         
   for :math:`|x| < 1`.  The range is extended to :math:`x < -1` by the
   Pfaff transformations [DLMF 15.8.1 and 15.8.2], which map the argument
   into :math:`(1/2, 1)`.  If the arguments :math:`(a,b,c,x)` are too close to a singularity then
   the function can return the error code :macro:`GSL_EMAXITER` when the
   series approximation converges too slowly.  This occurs in the region of
   :math:`x = 1`, :math:`c - a - b = m` for integer m.

   The contiguous cases :math:`c = a - 1` and :math:`c = b - 1` are
   evaluated from the closed forms

   .. math:: {}_2F_1(a,b,a-1,x) = {a - 1 + (b + 1 - a)x \over (a-1)(1-x)^{b+1}},

   .. math:: {}_2F_1(a,b,b-1,x) = {b - 1 + (a + 1 - b)x \over (b-1)(1-x)^{a+1}}.

.. exceptions:

.. function:: double gsl_sf_hyperg_2F1_conj (double aR, double aI, double c, double x)
              int gsl_sf_hyperg_2F1_conj_e (double aR, double aI, double c, double x, gsl_sf_result * result)

   These routines compute the Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, aR - i aI, c, x)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, aR - i aI, c, x)

   with complex parameters for :math:`|x| < 1`, extended to
   :math:`x < -1` by the same transformation.

.. function:: double gsl_sf_hyperg_2F1_renorm (double a, double b, double c, double x)
              int gsl_sf_hyperg_2F1_renorm_e (double a, double b, double c, double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a,b,c,x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a,b,c,x) / \Gamma(c)

   for :math:`|x| < 1`, extended to :math:`x < -1` as for
   :func:`gsl_sf_hyperg_2F1`.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F1_conj_renorm (double aR, double aI, double c, double x)
              int gsl_sf_hyperg_2F1_conj_renorm_e (double aR, double aI, double c, double x, gsl_sf_result * result)

   These routines compute the renormalized Gauss hypergeometric function

   .. only:: not texinfo

      .. math:: {}_2F_1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   .. only:: texinfo

      .. math:: 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c)

   for :math:`|x| < 1`, extended to :math:`x < -1` by the same
   transformation.
.. exceptions:

.. function:: double gsl_sf_hyperg_2F0 (double a, double b, double x)
              int gsl_sf_hyperg_2F0_e (double a, double b, double x, gsl_sf_result * result)

   These routines compute the hypergeometric function
   
   .. only:: not texinfo
      
      .. math:: {}_2F_0(a,b,x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x)

   The series representation is a divergent hypergeometric series.
   However, for :math:`x < 0` we have 

   .. only:: not texinfo

      .. math:: {}_2F_0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

   .. only:: texinfo

      .. math:: 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)

   For :math:`x = 0`, or for :math:`a = 0` or :math:`b = 0` with any
   :math:`x`, the value is :math:`1`.

.. exceptions: GSL_EDOM
