
**theta.h** -- Theta constants in genus 2
=========================================

This modules covers the computation of theta constants on the Siegel
upper half space in genus 2 only. This culminates in the computation
of theta constants in uniform quasi-linear time on the Siegel
fundamental domain. Some of the basic functionality (theta
characteristics, duplication formula) works for any *g* provided that
*2^{2g}* fits in an **ulong**.

**Theoretical warning**: it is not proven at this point that the
 Newton iterations used to compute theta constants on a compact set of
 the fundamental domain (**theta2_newton**) converge correctly. The
 correctness of many functions from now on relies on the assumption
 that **theta2_newton** produces correct output. This assumption holds
 in practice and will not be further mentioned.

Types and constants
-------------------

Theta characteristics are encoded as **ulongs** of bit length *2g*;
the leftmost (resp. rightmost) *g* bits encode a vector *a* (resp. b)
in *{0,1}^g*. The usual theta characteristic is then *(a/2, b/2)*.

We use the **acb_ptr** type to encode tuples of complex
numbers. Vectors of theta constants will have length *2^{2g}*, and are
indexed by their characteristics. These vectors will therefore contain
many zeroes corresponding to odd theta characteristics. A similar
numbering is used for vectors of complex numbers in Borchardt means.

Several constants are used for the tuning of Newton iterations:

::
   THETA_NEWTON_MINPREC

Newton iterations will not be performed below this precision.

::
   THETA_NEWTON_LOSS

Upper bound on the number of bits lost in each Newton iteration. This
approximately an upper bound on *log |df|*, where *f* denotes the
invertible function used in the Newton iterations.

::
   THETA_NEWTON_DERIVATIVE_OFFSET

Minor tweak to the small parameter *epsilon* used in finite
differences computations.

::
   THETA_NEWTON_BASEPREC

Approximate precision where theta constants will be computed using the
naive algorithm, before starting Newton iterations.

::
   THETA_NEWTON_Y1
   THETA_NEWTON_Y2MAX

Parameters in the function **theta2_unif** governing the choice
between the Newton and naive algorithm, and the number of duplications
to use.

::
   THETA_NEWTON_TOL_EXP

Tolerance parameter in the choice of Newton vs. naive algorithm.   
   
::
   THETA_DER_LOSS

Upper bound on the bits of precision lost when computing derivatives
of theta constants via finite differences in the Newton
algorithm. This is related to *THETA_NEWTON_LOSS*.
   
::
   BORCHARDT_VERBOSE
   THETA_VERBOSE

Set these constants to a positive number and recompile the library if
you wish messages to appear in **stdout**.


Borchardt means
---------------

We always assume *g=2*.

::
   int acb_sqrt_goodpos(acb_t r, const acb_t z, slong prec);

Set *r* to the principal square root of *z* as defined by
**acb_sqrt**, and returns whether the real part of *r* is certainly
positive.

::
   void borchardt_sqrt(acb_t r, const acb_t z, slong prec);

Set *r* to a square root of *z*. This uses **acb_sqrt**, and some
rescaling in the case where *z* intersects the negative real axis:
otherwise we would incur dramatic precision losses.

::
   void borchardt_root_ui(acb_t r, const acb_t z, ulong e, slong prec)

Set *r* to an *e* th root of *z*, using rescaling when *z* intersects
the negative axis as in **borchardt_sqrt**.

::
   int borchardt_step(acb_ptr b, acb_srcptr a, slong prec)

Set *b* to the result of a Borchardt step starting from *a* by taking
square roots in good position. If such a choice cannot be made, return 0.

::
   void borchardt_mean_m0(arb_t m0, acb_srcptr a, slong prec);

Set *m0* to a lower bound for the distance between the origin and an
affine half-plane containing all the entries of *a*.

::
   void borchardt_mean_M0(arb_t M0, acb_srcptr a, slong prec);

Set *M0* to an upper bound on the absolute values of all entries of
*a*.

::
   void borchardt_mean_Delta0(arb_t Delta0, acb_srcptr a, slong prec);

Set *Delta0* to an upper bound on the sum of all distances between
pairs of elements in *a*.

::
   int borchardt_mean_nb_steps_before_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

Compute an upper bound on the number of Borchardt steps in good
position starting from *a*, after which we are guaranteed to reach the
quadratic convergence area. This constant will likely be too large
compared to reality. Return 0 upon failure.

::
   void borchardt_mean_nb_steps_after_quad_conv(fmpz_t nb, acb_srcptr a, slong prec);

Compute an upper bound on the number of Borchardt steps in good
position starting from *a* that are necessary to obtain the Borchardt
mean with precision *prec*, assuming we are in the quadratic
convergence area.

::
   int borchardt_mean_quad_conv_is_reached(acb_srcptr a, slong prec);

Return whether we are certain that *a* lies in the quadratic
convergence area.

::
   int borchardt_mean(acb_t r, acb_srcptr a, slong prec);

Attempt to set *r* to the Borchardt mean of the vector *a* such that
all the Borchardt steps are in good position. Return 0 upon failure.

::
   void borchardt_excl_half_planes(arf_struct* b, const acb_t z, slong prec)

Compute the arguments of all complex numbers making an angle of more
than *pi/2* in absolute value with *z*. The result is a reunion of
intervals encoded in *b*.

::
   int borchardt_mean_invalid(acb_srcptr a, slong prec);

Return whether it is certain that no Borchardt step in good position
can be taken from *a*.


Theta characteristics
---------------------

These functions cover all genera *g* such that *2^{2g}* fits in an
**ulong**.

::
   ulong theta_char_get_a(ulong ch, slong g);
   ulong theta_char_get_b(ulong ch, slong g);
   ulong theta_char_set_ab(ulong a, ulong b, slong g);

Conversions between theta characteristics and vectors *a,b* as
described in the "Types and constants" section above.

::
   int theta_char_dot_product(ulong a, ulong b, slong g);

Compute the dot product between *a* and *b*.

::
   int theta_char_is_even(ulong ch, slong g);

Return whether *ch* is an even theta characteristic.

::
   slong theta_char_get_label_g2(ulong ch);
   ulong theta_char_set_label_g2(slong label);

Conversions between theta characteristics and Dupont labels when *g=2*.

::
   int theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g);
   int theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g);

Detection of Goepel quadruples and syzygous triples of theta characteristics.   


Theta constants
---------------

From now on, we fix *g=2*, and write an element *tau* of the Siegel
upper half space as *[tau1, tau3; tau3, tau2]*.

::
   void theta_duplication(acb_ptr th2_2tau, acb_srcptr th_tau, slong prec);

Apply the theta duplication formula to the vector *th_tau* (containing
at least the values of the 4 fundamental theta constants) and store
the result in *th2_2tau* (large enough for 16 complex values).

::
   int theta_0123_naive_B(fmpz_t B, const acb_mat_t tau, slong prec);

Compute an upper bound on the size of the square box necessary to
compute the fundamental theta constants at *tau* to precision *prec*
using the naive algorithm. Return 0 upon failure.

::
   int theta_0123_naive(acb_ptr th, const acb_mat_t tau, slong prec);

Attempt to evaluate the four fundamental theta constants at *tau*
using the naive algorithm. Return 0 upon failure. This can happen for
instance if *tau* overlaps the boundary of the Siegel upper half
space.

::
   int theta2_naive(acb_ptr th, const acb_mat_t tau, slong prec);

Compute the squares of all 16 theta constants at *tau*, using the
naive algorithm to compute fundamental theta constants at *tau/2* and
the duplication formula.

::
   int theta2_inverse(acb_mat_t tau, acb_srcptr th, slong prec);

Assuming *th* contains the values of squared theta constants at a
certain matrix *tau* in the Siegel fundamental domain, try to compute
such a matrix. Return 0 upon failure. This can happen if the input is
incorrect, or imprecise enough that we cannot decide if the relevant
Borchardt means are in good position.

::
   int theta2_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th, slong prec);

Same as **theta2_inverse**, but this time we leave *tau3^2* in the
entries outside the diagonal. This prevents dramatic precision losses
in Newton iterations when *tau3* is close to zero.

::
   int theta2_invalid(acb_srcptr th2, slong prec);

Return whether it is certain that *th2* does not contain the values of
squared theta constants at some matrix in the fundamental domain.

::
   int theta_0123half_diff_naive(acb_mat_t dth, const acb_mat_t tau, slong prec);
   int theta_0123half_inverse(acb_mat_t tau, acb_srcptr th_half, slong prec);
   int theta_0123half_inverse_no_sqrt(acb_mat_t tau, acb_srcptr th_half, slong prec);
   int theta_0123half_inverse_diff(acb_mat_t dtau, const acb_mat_t tau, acb_srcptr th_half,
				slong prec);

Subroutines for use inside the Newton iterations.

::
   int theta2_newton_step(acb_ptr th_half, const acb_mat_t tau, acb_srcptr th_half_approx,
		       slong prec);

Given a matrix *tau* in a certain compact set inside the fundamental
domain, and approximate values of the quotients
*theta_i(tau/2)/theta_0(tau/2)* at precision roughly *prec/2*, compute
the values of these quotients to precision roughly *prec* using a
Newton iteration.

::
   slong theta2_newton_start_prec(slong prec);

Compute the precision at which we start doing Newton iterations in
order to eventually reach *prec* with little overhead.

::
   int theta2_newton(acb_ptr th2, const acb_mat_t tau, slong prec);

Given a matrix *tau* in a certain compact set inside the fundamental
domain, attempt to compute all squared theta constants at *tau* up to
a common scalar factor using Newton iterations. Return 0 upon failure;
this can happen if the input is incorrect or too imprecise. **The
correctness of this function is heuristic**: cf. theoretical warning
at the beginning of this file.

::
   slong theta_newton_k2(acb_mat_t w, const acb_mat_t z, slong prec);
   slong theta_newton_k1(acb_mat_t w, const acb_mat_t z, slong prec);

Compute the number of duplication steps for use in **theta2_unif**.   

::
   int theta_use_naive(const acb_mat_t tau, slong prec);

Decide whether the matrix *tau* (assumed to lie in the fundamental
domain) is close enough to the cusp that the naive algorithm should be
used to compute theta constants at *tau*.

::
   int theta_use_newton(const acb_mat_t tau, slong prec);

Decide whether the matrix *tau* (assumed to lie in the fundamental
domain) is far enough from the cusp that the Newton algorithm should
be used to compute theta constants at *tau*.

::
   ulong theta_transform_image_char(fmpz_t epsilon, ulong ch, const fmpz_mat_t eta);

Compute the action of a symplectic matrix *eta* on the theta constant
of characteristic *ch*.

::
   void theta_transform_matrix(fmpz_mat_t res, const fmpz_mat_t eta);

Package the information provided by **theta_transform_image_char** in
a *16 x 3* matrix in the manner of Dupont's thesis.

::
   void theta_transform(acb_ptr th_eta, const fmpz_mat_t eta, acb_srcptr th, slong prec);

Assuming the vector *th* contains the values of theta constants at a
matrix *tau* in the Siegel upper half space up to a common scalar
factor, set *th_eta* to the values of theta constants at *eta tau* up
to a common scalar factor.

::
   void theta2_transform(acb_ptr th2_eta, const fmpz_mat_t eta, acb_srcptr th2, slong prec);

Same as **theta_transform** for squares of theta constants.

::
   int theta2_unif(acb_ptr th2, const acb_mat_t tau, slong prec);

Given a matrix *tau* in the fundamental domain, attempt to compute its
squared theta constants up to a common scalar factor. This functions
works in uniform quasi-linear time in *prec*, but can fail if the
input is too imprecise.

::
   int theta2_renormalize(acb_ptr th2, acb_srcptr th2_proj, slong prec);

Given the values of squared theta constants *th2_proj* at a matrix
*tau* lying in the fundamental domain, compute the actual values of
these theta constants and store them in *th2*. Return 0 on failure;
this can happen if the input is invalid or too imprecise.

::
   void theta2_randtest(acb_ptr theta2, flint_rand_t state, slong prec);

Generate a vector of 16 "random" complex numbers which contain the
values of squared theta constants up to some scalar factor at some
point in the Siegel upper half plane.


Derivatives of theta constants
------------------------------

Derivatives of theta constants with respect to the entries *tau1,
tau2, tau3* of a matrix *tau* are encoded as the entries of a *16 x 3*
complex matrix. Nb: in the cases where theta constants are computed up
to a common scalar factor, these derivatives are not uniquely defined.

::
   void theta_der_set_pert(arb_t eps, slong prec);

Set *eps* to a small perturbation for the computation of finite
differences.

::
   int theta_der_set_error(mag_t error, const acb_mat_t tau, slong prec);

Set *mag* to an upper bound on the error on derivatives of theta
constants, assuming we computed them from finite differences.

::
   int theta2_der_naive(acb_ptr th2_tau, acb_mat_t dth2_tau,
		      const acb_mat_t tau, slong prec);

Compute derivatives of squares of theta constants at *tau* to
precision roughly *prec/2*, using the naive algorithm and finite
differences.

::
   int theta_0123_der_naive(acb_ptr th, acb_mat_t dth,
			 const acb_mat_t tau, slong prec);

Same as **theta2_der_naive** for the values of the four fundamental
theta constants.

::
   void theta_der_duplication(acb_ptr th2_2tau, acb_mat_t dth2_2tau,
			   acb_srcptr th_tau, const acb_mat_t dth_tau,
			   slong prec);

Apply the formal differentiation of the theta duplication formula.

::
   int theta2_der_newton_step(acb_ptr th_half, acb_mat_t dth_approx,
			   const acb_mat_t tau, acb_srcptr th_half_approx,
			   slong prec);

Same as **theta2_newton_step**, but derivatives of theta constants at
*tau* at precision roughly *prec/2* are also part of the output.

::
   int theta2_der_newton(acb_ptr th2, acb_mat_t dth2, const acb_mat_t tau, slong prec);

Same as **theta2_newton**, but derivatives of theta constants at *tau*
at precision roughly *prec/2* are also part of the output.
