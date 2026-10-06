.. index::
   single: coupling coefficients
   single: 3-j symbols
   single: 6-j symbols
   single: 9-j symbols
   single: Wigner coefficients
   single: Racah coefficients

The Wigner 3-j, 6-j and 9-j symbols give the coupling coefficients for
combined angular momentum vectors.  Since the arguments of the standard
coupling coefficient functions are integer or half-integer, the
arguments of the following functions are, by convention, integers equal
to twice the actual spin value.  For information on the 3-j coefficients
see Abramowitz & Stegun, Section 27.9.  The functions described in this
section are declared in the header file :file:`gsl_sf_coupling.h`.

3-j Symbols
-----------

.. function:: double gsl_sf_coupling_3j (int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc)
              int gsl_sf_coupling_3j_e (int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc, gsl_sf_result * result)

   These routines compute the Wigner 3-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left(
         \begin{array}{ccc}
           ja & jb & jc \\
           ma & mb & mc
         \end{array}
         \right)

   .. only:: texinfo

      | ( ja jb jc )
      | ( ma mb mc )

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM

6-j Symbols
-----------

.. function:: double gsl_sf_coupling_6j (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf)
              int gsl_sf_coupling_6j_e (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, gsl_sf_result * result) 

   These routines compute the Wigner 6-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left\{
         \begin{array}{ccc}
           ja & jb & jc \\
           jd & je & jf
         \end{array}
         \right\}

   .. only:: texinfo

      | { ja jb jc }
      | { jd je jf }

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM

.. function:: double gsl_sf_coupling_6j_INCORRECT (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf)
              int gsl_sf_coupling_6j_INCORRECT_e (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, gsl_sf_result * result)

   .. deprecated:: 1.3

   These routines are deprecated and will be removed in a future release.
   They compute a permuted (incorrect) Wigner 6-j symbol and are retained
   only for backwards compatibility.  Use :func:`gsl_sf_coupling_6j` or
   :func:`gsl_sf_coupling_6j_e` instead.

9-j Symbols
-----------

.. function:: double gsl_sf_coupling_9j (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
              int gsl_sf_coupling_9j_e (int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji, gsl_sf_result * result) 

   These routines compute the Wigner 9-j coefficient, 

   .. only:: not texinfo

      .. math::

         \left\{
         \begin{array}{ccc}
           ja & jb & jc \\
           jd & je & jf \\
           jg & jh & ji
         \end{array}
         \right\}

   .. only:: texinfo
   
      | { ja jb jc }
      | { jd je jf }
      | { jg jh ji }

   where the arguments are given in half-integer units, :math:`ja` =
   :data:`two_ja`/2, :math:`ma` = :data:`two_ma`/2, etc.
.. Exceptional Return Values: GSL_EDOM

Wigner d-matrix
---------------

.. function:: double gsl_sf_wigner_drot (int two_j, int two_m1, int two_m2, double theta)
              int gsl_sf_wigner_drot_e (int two_j, int two_m1, int two_m2, double theta, gsl_sf_result * result)

   These routines compute the Wigner (small) d-matrix element

   .. math::

      d^{j}_{m_1 m_2}(\theta) = \langle j m_1 | e^{-i \theta J_y} | j m_2 \rangle

   where the angular momenta are given in half-integer units, :math:`j` =
   :data:`two_j`/2, :math:`m_1` = :data:`two_m1`/2, :math:`m_2` =
   :data:`two_m2`/2, and :data:`theta` is the rotation angle in radians.
   The value is computed from the finite sum of Zare, Eq. 3.57.  An
   error is returned if :math:`j < 0`, :math:`|m_1| > j`, :math:`|m_2| > j`
   or :math:`j + m_1` and :math:`j + m_2` are not integers.
.. Exceptional Return Values: GSL_EDOM
