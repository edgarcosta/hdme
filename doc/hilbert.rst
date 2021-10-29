
**hilbert.h** -- Hilbert and Humbert surfaces
=============================================

This modules covers computations on Hilbert and Humbert surfaces for
various small fundamental discriminants, and their link to the Siegel
moduli space.


Types and constants
-------------------

The constant *delta* always encodes the fundamental
discriminant. Elements in the ring of integers of the real quadratic
fields are encoded as **fmpz_poly_t** elements of degree at most 1;
the variable *w* in *aw+b* stands for the "fundamental generator", i.e.
*w = ((delta % 2) + sqrt(delta))/2*.

With our conventions, elements of the Hilbert modular group are
matrices *[a, b/sqrt(delta); c sqrt(delta), d]* where *a, b, c, d* are
quadratic integers. Such a matrix is encoded in an element of type
**fmpz_poly_mat_t** with *a, b, c, d* as entries in the above convention.

The type **acb_ptr** is used for periods in Hilbert space and
parameters for the various Hilbert or Humbert parametrizations.

::
   HILBERT_LLL_VERBOSE

Set this constants to a positive number and recompile the library if
you wish messages to appear in **stdout**.


Basic functionality for real quadratic fields
---------------------------------------------

::
   int hilbert_is_fundamental(slong delta);

Decide whether *delta* is a fundamental discriminant.

::
   int hilbert_is_totally_positive(const fmpz_poly_t x, slong delta);

Decide whether *x* defines a totally positive integer.

::
   int hilbert_splits(fmpz_poly_t beta, slong ell, slong delta);

Given a prime ell, run through small-trace quadratic integers in an
attempt to find a totally positive integer of norm *ell*. Upon
success, set *beta* to the result and return 1; otherwise leave *beta*
undefined and return 0.

::
   void hilbert_conjugate(fmpz_poly_t xbar, fmpz_poly_t x, slong delta);

Set *xbar* to the real conjugate of *x*.


Hilbert and Humbert parametrizations
------------------------------------

We require that *delta <= 100* for the Humbert parametrizations, and
*delta <= 20* for the Hilbert parametrizations. Gundlach invariants
are only implemented for *delta = 5*.

::
   void humbert_vars_set(char** vars, slong delta);

Store the names of the two variables used in the Humbert
parametrizations for discriminant *delta*.

::
   void humbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
		       slong delta, const fmpq_mpoly_ctx_t ctx);
   void gundlach_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
			slong delta, const fmpq_mpoly_ctx_t ctx);
   void hilbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
		       slong delta, const fmpq_mpoly_ctx_t ctx);

Set *pol* to the polynomial stored in the filed name *delta/name* in
the Humbert (resp. Gundlach, resp. Hilbert) folder.

::
   void humbert_AA1BB1B2(acb_ptr AA1BB1B2, acb_srcptr rs, slong delta,
		      slong prec);

Compute the quantities *A, A1, B, B1, B2* associated with the Humbert
parametrization of discriminant *delta* at the parameter *rs*; except
if *delta = 33*, in which case we directly compute the Igusa scalar
covariants.

::
   void humbert_cov_from_AA1BB1B2(acb_ptr I, acb_srcptr AA1BB1B2, slong prec);
   
Compute the Igusa scalar covariants associated with the given values
of *A, A1, B, B1, B2*; the same formula is used for all discriminants
*delta != 33*.

::
   void humbert_parametrize(acb_ptr I, acb_srcptr rs, slong delta, slong prec);
   void hilbert_parametrize(acb_ptr I, acb_srcptr rs, slong delta, slong prec);

Compute the value Igusa scalar covariants at the parameter *rs* in the
precomputed Humbert (resp. Hilbert) parametrization.
   

Gundlach invariants
-------------------

At present, these invariants are only implemented for *delta = 5*.

::
   void gundlach_from_igusa(acb_ptr g, acb_srcptr I, slong delta, slong prec);
   void gundlach_from_cov(acb_ptr g, acb_srcptr G, slong delta, slong prec);
   void cov_from_gundlach(acb_ptr G, acb_srcptr g, slong delta, slong prec);
   void gundlach_cov_from_igusa(acb_ptr G, acb_srcptr I, slong delta, slong prec);
   void igusa_from_gundlach(acb_ptr j, acb_srcptr g, slong delta, slong prec);

Conversions between Gundlach invariants *g1, g2*, Gundlach covariants
*G2, F6, F10*, Igusa scalar covariants and Igusa invariants.

::
   void gundlach_from_hilbert_param(fmpq* g, const fmpq* mn, slong delta);

Compute the Gundlach invariants associated with the parameter *mn*
given as rational numbers in the precomputed Hilbert parametrization
for *delta = 5*.
   

The Hilbert half space
----------------------

::
   void hilbert_halfspace_randtest(acb_ptr t, flint_rand_t state, slong prec);

Generate a "random" element *t* in the Hilbert upper half space whose
coefficients have approximately *prec* bits of precision. No claim
about the distribution is made.

::
   void hilbert_R(acb_mat_t R, slong delta, slong prec);

Compute a complex approximation of the standard matrix *R* used to
define the map from Hilbert space to Siegel space.

::
   void hilbert_map(acb_mat_t tau, acb_srcptr t, slong delta, slong prec);

Set *tau* to the image of *t* under the Hilbert embedding.

::
   int hilbert_linear_combination(fmpz* abcde, const acb_mat_t tau, slong delta, slong prec);

Given a matrix *tau* in Siegel space, attempt to use the LLL algorithm
to reconstruct integers *a, b, c, d, e* that are the coefficients of a
linear combination of the entries of *tau*, showing that the
corresponding abelian surface has real multiplication by
*Q(sqrt(delta))*. Return 0 upon failure.

::
   int hilbert_inverse(acb_ptr t, fmpz_mat_t eta, const acb_mat_t tau,
		    slong delta, slong prec);

Given a matrix *tau* in Siegel space, run the Birkenhake--Wilhelm
algorithm to compute a symplectic matrix *eta* such that *eta tau*
lies on the image of the Hilbert embedding, and set *t* to the
preimage of *eta tau*. Return 0 upon failure of
**hilbert_linear_combination**.

::
   void hilbert_sigma1(acb_t z, const fmpz_poly_t x, slong delta, slong prec);
   void hilbert_sigma2(acb_t z, const fmpz_poly_t x, slong delta, slong prec);

Set *z* to the image of *x* in the first (resp. second) complex
embedding of *Q(sqrt(delta))*.

::
   void hilbert_star(acb_t z, const fmpz_poly_mat_t m, acb_srcptr t,
		  slong delta, slong prec);

Set *z* to *(sigma1(c)t1 + sigma1(d), sigma2(c)t2 + sigma2(d))* where
*c, d* are the lower entries of *m* and *t = (t1, t2)*.

::
   void hilbert_transform(acb_ptr z, const fmpz_poly_mat_t m, acb_srcptr t,
		       slong delta, slong prec);

Set *z* to the image of *t* under the Hilbert modular transformation *m*.

::
   void hilbert_scalar_mul(acb_ptr z, const fmpz_poly_t x, acb_srcptr t,
		        slong delta, slong prec);

Set *z* to *(sigma1(x) t1, sigma2(x) t2)*.			

::
   void hilbert_scalar_div(acb_ptr z, const fmpz_poly_t x, acb_srcptr t,
			slong delta, slong prec);

Set *z* to *(t1/sigma1(x), t2/sigma2(x))*.

::
   void hilbert_transform_randtest(fmpz_poly_mat_t m, flint_rand_t state, slong bits);

Set *m* to a "random" matrix in the Hilbert modular group whose
coefficients have size approximately *bits*. No claim about the
distribution is made.

