**siegel.h** -- The Siegel upper half space
===========================================

This module covers computations with complex matrices in the Siegel
upper half space for genus 2, culminating with reduction to the
classical fundamental domain for the action of the modular group
*Sp_4(ZZ)*. Some of the basic functions support arbitrary genus.

Types and constants
-------------------

We use the type **acb_mat_t** for elements of the Siegel half space,
and the type **fmpz_mat_t** for symplectic matrices.

Additional functions for real and complex matrices
--------------------------------------------------

::
   void acb_mat_get_real(arb_mat_t re, const acb_mat_t z);
   void acb_mat_get_imag(arb_mat_t im, const acb_mat_t z);

Set *re* and *im* to the real and imaginary part of *z* respectively.

::
   void acb_mat_set_arb_arb(acb_mat_t z, const arb_mat_t re, const arb_mat_t im);

Set *z* to the complex matrix with real and imaginary parts *re*, *im*.

::
   void arb_mat_randtest_precise(arb_mat_t r, flint_rand_t state, slong prec,
			      slong mag_bits);

Set the coefficients of *r* to "random" real numbers of precision
*prec* and magnitude given by *mag_bits* approximately, as described
in the **arb_randtest_precise** function from Arb.

::
   void arb_mat_randtest_sym_precise(arb_mat_t r, flint_rand_t state, slong prec,
				  slong mag_bits);

Same as **arb_mat_randtest_precise**, with the additional property
that *r* is symmetric.

::
   void arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec);

Set *r* to a "random" symmetric and positive definite matrix of
precision approximately prec. It is guaranteed that the box *r*
contains a symmetric positive definite matrix.

::
   int arb_mat_is_nonsymmetric(const arb_mat_t m);

Return whether *m* is certainly not symmetric.

::
   void arb_mat_congr_fmpz_mat(arb_mat_t r, const fmpz_mat_t u, const arb_mat_t m,
			    long prec);

Set *r* to *u m u^t*.

::
   int arb_mat_is_minkowski_reduced(const arb_mat_t r, const arb_t tol, slong prec);

Return whether it is certain that the matrix *r* is approximately
Minkowski-reduced with tolerance parameter *tol*. See `[K21]`_ for the
precise definition. This functions assumes *g <= 2*.

::
   int arb_mat_not_minkowski_reduced(const arb_mat_t r, slong prec);

Return whether it is certain that *r* is not Minkowski-reduced. This
function assumes *g <= 2*.

::
   int arb_mat_minkowski_reduce(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m,
			 const arb_t tol, slong prec);

Attempt to compute an approximate Minkowski reduction of *m* with
tolerance parameter *tol*; set *r* to the result and *u* to a matrix
in *SL_2(ZZ)* such that *r = u m u^t*. For best results, *tol* should
be order of magnitudes larger than *2^{-prec}*. Return 1 upon success
and 0 upon failure. This function assumes *g <= 2*.

::
   void arb_mat_lambda(arb_t lambda, const arb_mat_t m, slong prec);
   
Compute a lower bound for the eigenvalues of the imaginary part *im*
of *m*, assuming that *im* is positive definite. This function assumes
*g <= 2*.


Symplectic matrices
-------------------

::
   slong fmpz_mat_half_dim(const fmpz_mat_t m);

Return *g* such that *m* is a square matrix of size *2g x 2g*. Abort
if *m* is not a square matrix of even dimension.

::
   void fmpz_mat_direct_inv(fmpz_mat_t minv, const fmpz_mat_t m);

Set *minv* to the inverse of *m*. Abort if *m* is not invertible.

::
   void fmpz_mat_get_a(fmpz_mat_t a, const fmpz_mat_t m);
   void fmpz_mat_get_b(fmpz_mat_t a, const fmpz_mat_t m);
   void fmpz_mat_get_c(fmpz_mat_t a, const fmpz_mat_t m);
   void fmpz_mat_get_d(fmpz_mat_t a, const fmpz_mat_t m);

Set *a* to the upper-left (res. upper-right, lower-left, lower-right)
*g x g* block of *m*.

::
   void fmpz_mat_set_abcd(fmpz_mat_t m,
		       const fmpz_mat_t a, const fmpz_mat_t b,
		       const fmpz_mat_t c, const fmpz_mat_t d);

Set *m* to the matrix with *g x g* blocks *a, b, c, d*.		       

::
   void fmpz_mat_J(fmpz_mat_t m);

Set *m* to the symplectic matrix *J* with blocks *0, I_g, -I_g, 0*.

::
   int fmpz_mat_is_J(const fmpz_mat_t m);

Return whether *m = J*.

::
   int fmpz_mat_is_symplectic(const fmpz_mat_t m);

Return whether *m* is symplectic. Abort if *m* is not a square matrix
of even dimension.

::
   void fmpz_mat_diagonal_symplectic(fmpz_mat_t m, const fmpz_mat_t u);

Set *m* to the symplectic matrix with blocks *u, 0, 0, u^{-t}*. Abort
if *u* is not invertible.

::
   void fmpz_mat_randtest_triangular_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);

Set *m* to a "random" symplectic matrix with blocks *I_g, S, 0, I_g*,
where *S* is a "random" symmetric integer matrix with coefficients of
size approximately *bits*.

::
   void fmpz_mat_randtest_diagonal_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);

Set *m* to a "random" symplectic matrix with blocks *u, 0, 0, u^{-t}*
with coefficients of size approximately *bits*.

::
   void fmpz_mat_randtest_symplectic(fmpz_mat_t m, flint_rand_t state, slong bits);

Set *m* to a "random" symplectic matrix with coefficients of size
approximately *bits*, using a combination of the functions above. No
claim on the distribution is made.


The Siegel upper half space
---------------------------

::
   void siegel_halfspace_randtest(acb_mat_t z, flint_rand_t state, slong prec);

Set *z* to a "random" element of the Siegel upper half space whose
coefficients have precision approximately *prec* and magnitude
*O(1)*. No claim of the distribution is made.

::
   void siegel_star(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec);

Set *w* to *cz+d*, where *c, d* are the lower blocks of *m*.

::
   int siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec);

Set *w* to *(az+b)(cz+d)^{-1}*, where *a,b,c,d* are the blocks of
*m*. Return *0* upon failure to compute the inverse.

::
   int siegel_is_real_reduced(const acb_mat_t z, const arb_t tol, slong prec);

Return whether it is certain that the real part of *z* is
approximately reduced with tolerance parameter *tol*. See `[K21]`_ for
the precise definition.

::
   int siegel_not_real_reduced(const acb_mat_t z, slong prec);

Return whether it is certain that the real part of *z* is not reduced.   

::
   int siegel_reduce_real(acb_mat_t w, fmpz_mat_t u, const acb_mat_t z,
		       const arb_t tol, slong prec);

Attempt to reduce the real part of *z* with tolerance parameter *tol*;
set *w* to the result and *u* to a symplectic matrix such that *w =
uz*. For best results, *tol* should be order of magnitudes larger than
*2^{-prec}*. Return 1 upon success and 0 upon failure.

::
   slong siegel_nb_test_matrices(slong g);

Return the number of symplectic matrices defining the boundary of the
fundamental domain in the Siegel half space of genus *g*. This
functions assumes *g <= 2*.

::
   void siegel_test_matrix(fmpz_mat_t u, slong j);

Set *u* to the test matrix number *j*. This function assumes *g <= 2*.

::
   int siegel_fundamental_domain(acb_mat_t w, fmpz_mat_t m,
			      const acb_mat_t z, const arb_t tol, slong prec);

Attempt to reduce *z* to a neighborhood of the Siegel fundamental
domain, with tolerance parameter *tol*; set *w* to the result and *m*
to a symplectic matrix such that *w = mz*. For best results, *tol*
should be order of magnitudes larger than *2^{-prec}*. Return 1 upon
success and 0 upon failure. This function assumes *g <= 2*. See
`[K21]`_ for the precise definition of this neighborhood.

::
   int siegel_is_in_fundamental_domain(const acb_mat_t z, const arb_t tol, slong prec);

Return whether it is certain that *z* lies in the neighborhood of the
fundamental domain specified by the tolerance parameter *tol*. This
function assumes *g <= 2*.

::
   int siegel_not_in_fundamental_domain(const acb_mat_t z, slong prec);

Return whether it is certain that *z* does not lie in the Siegel
fundamental domain. This function assumes *g <= 2*.

::
   int siegel_is_weakly_reduced(const acb_mat_t z, const arb_t tol, slong prec);

Return whether it is certain that *z* lies in the neighborhood of the
enlarged domain *F_2'*, specified by the tolerance parameter
*tol*. See `[K21]`_ for the precise definition. This function assumes
*g <= 2*.

::
   void siegel_fundamental_domain_randtest(acb_mat_t z, flint_rand_t state, slong prec);

Set *z* to a "random" matrix in the Siegel fundamental domain whose
coefficients have precision approximately *prec* and magnitude *O(1)*.

.. _[K21]: https://arxiv.org/abs/2010.10094
