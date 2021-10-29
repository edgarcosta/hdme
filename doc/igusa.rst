
**igusa.h** -- Igusa invariants and Mestre's algorithm
======================================================

This modules covers the conversions between theta constants in genus
2, Igusa scalar covariants, Igusa invariants. Conversions to the
coefficients of an associated genus 2 hyperelliptic curve via Mestre's
algorithm and Thomae's formulae are also part of this module.

We use Streng's version of Igusa invariants.

Types and constants
-------------------

The type **acb_ptr** is systematically used to record vectors of
invariants, coefficients, roots, etc.

::
   THOMAE_LOWPREC

Minimum precision used to rule out bad sign choices in Thomae's
formulae.

::
   THOMAE_MULPREC

Gap between the precisions used at successive step to refine the
possible good sign choices in Thomae's formulae.
   
::
   THOMAE_VERBOSE
   
Set this constants to a positive number and recompile the library if
you wish messages to appear in **stdout**.


Igusa covariants from theta constants
-------------------------------------

::
   void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec);
   void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec);
   void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec);
   void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec);
   void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec);

Compute the values of the usual Siegel modular forms *h4, h6, h10,
h12, h16* given the values of squares of theta constants stored in
*theta2*.

::
   void igusa_h(acb_ptr h, acb_srcptr theta2, slong prec);

Set the four first entries of *h* to *h4, h6, h10, h12*.

::
   void cov_from_h(acb_ptr I, acb_srcptr h, slong prec);

Compute the values of the Igusa scalar covariants *I2, I4, I6', I10*
associated with the given values of *h4, h6, h10, h12*.

::
   void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec);

Compute the values of the Igusa scalar covariants associated with the
given values of squared theta constants.

::
   void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal, slong prec);

Set *I* to the vector *S* rescaled by *scal^k* for *k = 2, 4, 6, 10*.

::
   void igusa_from_theta2(acb_ptr j, acb_srcptr theta2, slong prec);

Compute the values of Igusa invariants (Streng's version) from the
given values of squares of theta constants.

::
   int igusa_from_tau(acb_ptr j, const acb_mat_t tau, slong prec);

Compute the values of Igusa invariants at the given matrix *tau*
assuming it lies in the Siegel fundamental domain. Return 0 when the
call to **theta2_unif** fails.

::
   int igusa_is_defined(acb_srcptr j);

Return whether the three Igusa invariants are all finite complex
numbers.


Igusa covariants from curve coefficients
----------------------------------------

Genus 2 hyperelliptic curves are encoded as degree 6 polynomials.

::
   void curve_coeffs(acb_ptr ai, const acb_poly_t crv);
   void curve_coeffs_fmpz(fmpz* ai, const fmpz_poly_t crv);

Store the coefficients of *crv* in the vector *ai*.

::
   void igusa_scalar_covariants(acb_ptr I, const acb_poly_t crv, slong prec);
   void igusa_scalar_covariants_fmpz(fmpz* I, const fmpz_poly_t crv);

Compute the values of Igusa scalar covariants *I2, I4, I6, I10* at the
given curve equation.

::
   void igusa_from_cov(acb_ptr j, acb_srcptr I, slong prec);
   void igusa_from_cov_fmpz(fmpq* j, const fmpz* I);

Compute the values of Igusa invariants (Streng's version) from the
given values of Igusa scalar covariants.
   
::
   void cov_from_igusa(acb_ptr I, acb_srcptr j, slong prec);

Compute possible values of Igusa scalar covariants giving rise to the
given Igusa invariants.

::
   void igusa_from_curve(acb_ptr j, const acb_poly_t crv, slong prec);
   void igusa_from_curve_fmpz(fmpq* j, const fmpz_poly_t crv);

Compute the Igusa invariants of the given genus 2 hyperelliptic curve.


Different covariants: I6, and Clebsch
-------------------------------------

::
   void igusa_I6(acb_t I6, acb_srcptr I, slong prec);

Compute the value of the scalar covariant *I6* from the given values
of *I2, I4, I6'*.

::
   void igusa_clebsch(acb_ptr ABCD, acb_srcptr I, slong prec);

Compute the Clebsch invariants *A,B,C,D* in Mestre's notation from the
given values of *I2, I4, I6', I10*.

::
   void igusa_R2(acb_t res, acb_srcptr I, slong prec);

Compute the value of the degree-30 covariant *R^2* from the given
values of *I2, I4, I6, I10*.


Mestre's algorithm
------------------

Since we work with inexact input, only the generic version of Mestre's
algorithm is implemented.

::
   int igusa_has_generic_automorphisms(acb_srcptr I, slong prec);

Return whether it is certain that any curve with the given Igusa
scalar covariants has only *{+/-1}* as automorphisms.

::
   void igusa_generic_randtest(acb_poly_t crv, acb_ptr I, flint_rand_t state, slong prec);

Generate a "random" genus 2 curve with generic automorphisms and
compute its Igusa covariants.

::
   void mestre_conic(acb_ptr conic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

Store the coefficients of (a rescaled version of) Mestre's conic in
the vector *conic*, given *A, B, C, D, U*, and *I10 = D'* in Mestre's
notation.

::
   void mestre_conic_randtest(acb_ptr conic, flint_rand_t state, slong prec);

Generate a "random" tuple of coefficients looking like the
coefficients of Mestre's conic.

::
   void mestre_cubic(acb_ptr cubic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

Store the coefficients of (a rescaled version of) Mestre's cubic in
the vector *cubic*.

::
   void mestre_subst_in_conic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr conic, slong prec);
   void mestre_subst_in_cubic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr cubic, slong prec);

Set *subst* to the evaluation of Mestre's conic (resp. cubic) at the
triple *(x1, x2, x3)*.

::
   int mestre_point_on_conic(acb_ptr pt, acb_srcptr conic, slong prec);

Set *pt* to a box containing a point on Mestre's conic. This function
can fail if the input is too imprecise.

::
   int mestre_point_is_outside_conic(acb_srcptr pt, acb_srcptr conic, slong prec);

Return whether it is certain that the given point lies outside
Mestre's conic.

::
   void mestre_eval_cubic(acb_t res, acb_srcptr pt, acb_srcptr cubic, slong prec);

Set *res* to the evaluation of Mestre's cubic at the given complex
point *pt*.

::
   void mestre_parametrize_conic(acb_poly_t x1, acb_poly_t x2, acb_poly_t x3,
			     acb_srcptr pt, acb_srcptr conic, slong prec);

Set *x1, x2, x3* to polynomials parametrizing Mestre's conic, given a
base point *pt* lying on it.

::
   int mestre(acb_poly_t crv, acb_srcptr I, slong prec);

Reconstruct a genus 2 curve *crv* realizing the given Igusa scalar
covariants. Return 0 upon failure, e.g. when the curve is not
guaranteed to have generic automorphisms or the input is too
imprecise.
   

Thomae's formulae
-----------------

When applying Thomae's formula, a choice of ordering of the complex
roots of *crv* and certain sign choices have to be made. We encode
both of these choices as **slongs**: the permutation number between 1
and 720, and a flag between 0 and 15 for the sign choices.

::
   int thomae_roots(acb_ptr roots, const acb_poly_t crv, slong prec);

Set *roots* to the vector of complex roots of the polynomial *crv*,
using **acb_poly_find_roots** as a subroutine. Return 0 if the roots
could not be properly isolated, e.g. if the curve is not guaranteed to
be nonsingular.

::
   void thomae_reorder(acb_ptr new_roots, acb_srcptr roots, slong perm);

Reorder the roots according to the permutation number *perm*.

::
   void thomae_rosenhain(acb_ptr ros, acb_srcptr roots, slong prec);

Compute the three Rosenhain invariants associated with the given
ordered set of roots.

::
   void thomae_theta4(acb_ptr th4, acb_srcptr ros, slong prec);

Compute all 16 fourth powers of theta constants corresponding to the
given values of Rosenhain invariants.

::
   void thomae_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros, slong signs, slong prec);

Compute all 16 squares of theta constants corresponding to the given
Rosenhain invariants and choice of signs.

::
   int thomae_discard(acb_srcptr th2, slong prec);

Return 1 if the given vector of theta constants certainly does not
come from the correct ordering of roots and choice of signs.

::
   int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec);

Return 0 if it is certain that the candidate period matrix *tau* is
not a matrix in the Siegel fundamental domain with the prescribed
values of Igusa scalar covariants up to rescaling.

::
   slong thomae_startprec(slong prec);

Compute the precision at which iterations of **thomae_correct_signs**
should start in order to finally reach *prec* with little overhead.

::
   int thomae_correct_signs(slong* perm, slong* signs, acb_srcptr roots,
			 acb_srcptr I, slong prec);

Attempt to compute the correct permutation and sign choices from the
given set of roots, in order to obtain a period matrix in the Siegel
fundamental domain realizing the given Igusa scalar covariants up to
rescaling. Return 0 upon failure, e.g. if the input is incorrect or
too imprecise.

::
   int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec);

Attempt to compute a period matrix in the Siegel fundamental domain
realizing the given Igusa scalar covariants up to rescaling, by
combining Mestre's algorithm with Thomae's formulae.

::
   int tau_theta2_from_igusa(acb_mat_t tau, acb_ptr th2, acb_srcptr I, slong prec);

Same as **tau_from_igusa**, but the values of renormalized squared
theta constants at *tau* are also part of the output.

#endif 
