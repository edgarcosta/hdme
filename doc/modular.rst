
**modular.h** -- Evaluation of modular equations in genus 2
===========================================================

This modules covers the evaluation of modular equations in genus 2 at
a given point over either the field of rational numbers or a prime
finite field. This combines functions from all other modules.

The types of modular equations currently covered are:

1. Siegel modular equations in Igusa invariants.

2. Symmetric Hilbert modular equations in Igusa invariants, evaluated
  at a given point in the precomputed Humbert parametrization for all
  fundamental discriminants *delta <= 100*.

3. Nonsymmetric Hilbert modular equations in Igusa invariants,
   evaluated at a given point in the precomputed Hilbert
   parametrization for all fundamental discriminants *delta <= 20*.

4. Symmetric Hilbert modular equations in Gundlach invariants for
   *delta = 5*.

5. Nonsymmetric Hilbert modular equations in Gundlach invariants for
   *delta = 5*, evaluated at a given point in the precomputed Hilbert
   parametrization.

Assuming correct output from **theta2_newton**
(resp. **theta2_newton** and **hilbert_linear_combination**), the
evaluation of modular equations of type 1. (resp. types 4., 5.)
returns a provably correct result. Otherwise, there is a heuristic
rational reconstruction step.

Types and constants
-------------------

Modular equations of types 1, 2, 3 are given by three polynomials;
modular equations of types 4, 5 are given by two polynomials only. In
both cases we store the numerators in a vector with entries of type
**fmpz_poly_struct**, and the denominator as a separate **fmpz_t**.

::
   MODEQ_VERBOSE

Set this constant to 0 and recompile the library if you do not wish
messages to appear in **stdout**.

::
   MODEQ_RED_TOL_BITS

Tolerance in reduction to the Siegel fundamental domain will be *1/2*
to this power.

::
   MODEQ_MAX_PREC

Maximal precision allowed in the computation of modular
equations. This bound prevents infinite computations if rational
coefficients cannot be recognized. For larger examples, it is
possible that the default value has to be increased.

::
   SIEGEL_START_PREC_MUL
   SIEGEL_START_PREC_ADD
   SIEGEL_MUL_PREC

Constants governing the choice of successive precisions in the
algorithm for Siegel modular equations.

::
   HILBERT_START_PREC_MUL
   HILBERT_START_PREC_ADD
   HILBERT_MUL_PREC

Constants governing the choice of successive precisions in the
algorithms for the various types of Hilbert modular equations. The
default values are particularly adapted to the computation of
symmetric Hilbert modular equations in Gundlach invariants for *delta
= 5*. Expect a large number of precision increases in other cases.


All types of modular equations
------------------------------

::
   void product_tree_1(acb_poly_t P, acb_srcptr xi, acb_srcptr yi, slong d, slong prec);

Set *P* to the polynomial *\prod_{i=1}^d (x_i X + y_i)* using a
product tree.

::
   void product_tree_2(acb_poly_t Q, acb_srcptr xi, acb_srcptr yi, acb_srcptr zi,
		    slong d, slong prec);

Set *Q* to the polynomial *\sum_i z_i \prod_{j\neq i} (x_j X + y_j)*
using another product tree.

::
   int modeq_round_coeff(fmpz_t c, const acb_t x);

If the ball defined by *x* contains a unique integer, set *c* to this
value and return 1; otherwise leave *c* undefined and return 0.

::
   int modeq_round_poly(fmpz_poly_t pol, arf_t max_radius,
		     const acb_poly_t pol_acb, slong degree);

If the ball defined by *pol_acb* contains a unique polynomial with
integer coefficients, set *pol* to this value and return 1; otherwise
return 0. In any case *max_radius* is set to the maximum radius found
among all the coefficients of *pol_acb*.

::
   void modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong nb, slong prec);

Given a vector *th2_vec* of length *16 nb* containing *nb* 16-tuples
of squared theta constants, fill *I_vec* with *nb* 4-tupes of
associated Igusa scalar covariants.

::
   int modeq_round(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* num_vec_acb,
		const acb_t den_acb, slong degree, slong nb);

Attempt to recognize *nb* polynomials with integer coefficients from
*num_vec_acb*, and an integer from *den_acb*. In case of success, set
*num_vec* and *den* to the results and return 1; otherwise return 0.

::
   int modeq_rational_coeff(fmpq_t c, fmpz_t den, const acb_t x,
			 const fmpz_t probable_den, slong prec);

Attempt to recognize a rational number from the given complex value
*x* with high degree of certainty using a continued fraction
algorithm. The nonzero integer *probable_den* is used as a candidate
denominator to speed up the computation. If successful, set *c* to the
result and *den* to the lcm of the denominator of *c* and
*probable_den*; otherwise leave these undefined and return 0.

::
   int modeq_rational_poly(fmpq_poly_t pol, const acb_poly_t pol_acb,
			slong degree, slong prec);

Run **modeq_rational_coeff** on each coefficient of *pol_acb*. If a
polynomial with rational coefficients was recognized, set *pol* to the
result and return 1; otherwise return 0.

::
   int modeq_rational(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* pol_vec_acb,
		   slong degree, slong nb, slong prec);

Attempt to recognize *nb* polynomials with rational coefficients from
*pol_vec_acb*. If successful, set *den* to a common denominator of the
results and *num_vec* to the collection of the corresponding
numerators; otherwise return 0.

::
   void modeq_input_get_fmpz(fmpz_t den, fmpz* num, fmpq* j, slong len);

Set *den* to a common denominator of all entries of the vector *j* of
length *len*, and fill *num* with the corresponding numerators.

::
   slong modeq_height_fmpz(const fmpz* j, slong len);
   slong modeq_height_fmpq(fmpq* j, slong len);

Return an approximation of the logarithmic height of the entries of
*j*.

::
   void modeq_simplify(fmpz_poly_struct* num_vec, fmpz_t den, slong degree, slong nb);

Divide the *nb* entries of *num_vec* and the integer *den* by their
greatest common divisor in-place.

::
   void modeq_factor_Q(slong* nb_factors, fmpz_poly_struct* factors, slong* exps,
		    const fmpz_poly_t pol);

Factor the given polynomial *pol* by calling **fmpz_poly_factor**. Set
*nb_factor* to the number of factors, *factors* to the list of all
factors, and *exps* to the list of exponents.

::
   void modeq_roots_Q(slong* nb_roots, fmpq* roots, slong* mults,
		   const fmpz_poly_t pol);

Find the roots of the given polynomial *pol* by calling
**modeq_factor_Q**. Set *nb_roots* to the number of roots, *roots* to
the list of roots, and *mults* to the list of multiplicities.

::
   int modeq_isog_invariants_Q(fmpq* j, const fmpz_poly_struct* num_vec,
			    const fmpq_t root, slong nb);

Given the vector *num_vec* containing (numerators of) evaluated
modular equations and a root *root* of the first entry *num1* of
*num_vec* (encoding the first invariant of an abelian surface *B*
isogenous to the abelian surface *A* where we evaluated modular
equations), set *j* to the complete list of invariants of *B*. Return
0 on failure, i.e. if we hit the rare case where *num1* has a double
root at *j*.

::
   void modeq_factor_Fp(slong* nb_factors, fmpz_mod_poly_struct* factors, slong* exps,
		     const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);		     
   void modeq_roots_Fp(slong* nb_roots, fmpz* roots, slong* mults,
		    const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);
   int modeq_isog_invariants_Fp(fmpz* j, const fmpz_mod_poly_struct* pol_vec,
			     const fmpz_t root, slong nb,
			     const fmpz_mod_ctx_t ctx);

Same as above over the prime finite field encoded in the context *ctx*.			     

::
   void modeq_input_lift(fmpq* j, const fmpz* input, slong nb);

Set *j* to a copy of the vector *input* of length *nb*.

::
   int modeq_reduce(fmpz_mod_poly_struct* red_vec, const fmpz_poly_struct* num_vec,
		 const fmpz_t den, slong nb, const fmpz_mod_ctx_t ctx);

Given evaluated modular equations encoded in *num_vec* and *den* over
the rationals, reduce them to the finite field encoded in the context
*ctx* and store the result in *red_vec*. Return 0 upon failure,
i.e. *den* is zero in the finite field.


Siegel modular equations
------------------------

::
   slong siegel_modeq_height_fmpz(const fmpz* j);
   slong siegel_modeq_height_fmpq(fmpq* j);

Return an approximation of the logarithmic height of the entries of
*j*. These functions are customized to reflect the disparity between
partial degrees in modular equations of Siegel type.

::
   slong siegel_nb_cosets(slong ell);

Return the degree of the first evaluated Siegel modular equation,
i.e. *ell^3 + ell^2 + ell + 1*.

::
   void siegel_coset(fmpz_mat_t m, slong k, slong ell);

Set *m* to the *k* th element of the list of coset representatives for
the modular subgroup *Gamma^0(ell)* inside *Gamma(1)*, after
multiplication of the lower blocks by *ell* and transformation such
that the cusp at infinity is left stable.

::
   int siegel_modeq_theta2(acb_ptr th2_vec, acb_ptr stardets,
			const acb_mat_t tau, slong ell, slong prec);

Set *th2_vec* to the vector of normalized squared theta constants at
reductions to the fundamental domain of all period matrices which are
*ell*-isogenous to *tau*. The rescaling factors coming from the
reduction matrix (one for each isogenous period matrix) are stored in
the vector *stardets*. Return 0 upon failure, i.e. if one of the
computations of theta constants or one of the reductions fails.

::
   void siegel_modeq_exps(slong* e, slong* a, slong* b, slong ell);

Set *e, a, b* such that the scalar factor by which we multiply complex
evaluations of Siegel modular equations should be multiplied by
*I4(tau)^e / (I10(tau)^a I4(tau)^b)*.

::
   void siegel_modeq_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
			 slong ell, slong prec);

Compute a rescaling factor for complex evaluations of Siegel modular
equations such that we will be able to recognize integer coefficients.

::
   void siegel_modeq_num(acb_poly_struct* num_vec_acb,
		      acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);
   void siegel_modeq_den(acb_t den, acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);

Compute complex evaluations of Siegel modular equations using product
trees from the input vector *I_vec*, and multiply the results by *scal*.

::
   slong siegel_modeq_startprec_fmpz(const fmpz* j, slong ell);
   slong siegel_modeq_startprec_fmpq(fmpq* j, slong ell);

Choose a starting precision for the evaluation of Siegel modular
equations of level *ell* at the vector *j*. Usually, this starting
precision will be just large enough to conclude.

::
   slong siegel_modeq_nextprec(slong current_prec);

Choose the next precision level, in case the computation at precision
*current_prec* was not successful.

::
   void siegel_modeq_rescale(fmpz_t scal, fmpq* j, slong ell);

Compute another rescaling factor by which we should multiply Siegel
modular equation in order to recognize integers, in case the vector
*j* does not consist of integers to begin with.

::
   int siegel_modeq_eval_Q(fmpz_poly_struct* num_vec,
			fmpz_t den, fmpq* j, slong ell);
   int siegel_modeq_eval_Fp(fmpz_mod_poly_struct* pol_vec,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx);

Evaluate Siegel modular equations of level *ell* at the given vector
of Igusa invariants *j*. Return 0 upon failure, i.e. if the
computation did not succeed before reaching the maximal allowed
precision.


Hilbert modular equations
-------------------------

We always assume that the prime *ell* splits in the real quadratic
field *Q(sqrt(delta))* as the product of two principal ideals
generated by totally positive elements *beta* and its conjugate
*betabar*.

::
   slong hilbert_nb_cosets(slong ell, slong delta);

Return the degree of nonsymmetric Hilbert modular equations of level
*ell*, i.e. *ell+1*.

::
   void hilbert_coset(fmpz_poly_mat_t m, slong k, slong ell, slong delta);

Set *m* to the *k* th element in the list of coset representatives for
the subgroup *Gamma^0(beta)* inside the Hilbert modular group.

::
   int hilbert_modeq_theta2(acb_ptr th2_vec, acb_srcptr t,
			 const fmpz_poly_t beta, slong ell, slong delta, slong prec);
   int hilbert_modeq_theta2_star(acb_ptr th2_vec, acb_ptr stardets,
			      acb_srcptr t,
			      const fmpz_poly_t beta, slong ell, slong delta, slong prec);

Same as **siegel_modeq_theta2** for *beta*-isogenous periods in
Hilbert space instead of *ell*-isogenous periods in Siegel
space. Rescaling factors are only recorded in the **_star** version.

::
   void hilbert_modeq_igusa_C(acb_poly_struct* pol_vec,
			   acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
			   slong ell, slong delta, slong prec);
   void hilbert_modeq_nonsym_igusa_C(acb_poly_struct* pol_vec, acb_srcptr I_vec, slong ell,
				  slong delta, slong prec);

Compute complex evaluations of symmetric (resp. nonsymmetric) Hilbert
modular equations in Igusa invariants, given the relevant vectors of
Igusa scalar covariants at isogenous periods.

::
   void hilbert_modeq_gundlach_exps(slong* e, slong* a, slong* b, slong ell, slong delta);

Same as **siegel_modeq_exps** in the case of Hilbert modular equations
in Gundlach invariants for discriminant 5; this time the fraction
takes the form *G_2(tau)^e / ((F10(tau)^a G2(tau)^b))*.

::
   void hilbert_modeq_gundlach_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
				   slong ell, slong delta, slong prec);
   void hilbert_modeq_gundlach_num(acb_poly_struct* num_vec_acb,
				acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
				const acb_t scal,
				slong ell, slong delta, slong prec);
   void hilbert_modeq_gundlach_den(acb_t den, acb_srcptr I_vec_beta,
				acb_srcptr I_vec_betabar, const acb_t scal,
				slong ell, slong delta, slong prec);
   void hilbert_modeq_nonsym_gundlach_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
					  slong ell, slong delta, slong prec);
   void hilbert_modeq_nonsym_gundlach_num(acb_poly_struct* num_vec_acb,
				       acb_srcptr I_vec_beta,
				       const acb_t scal,
				       slong ell, slong delta, slong prec);
   void hilbert_modeq_nonsym_gundlach_den(acb_t den, acb_srcptr I_vec_beta,
				       const acb_t scal,
				       slong ell, slong delta, slong prec);

Analogues of **siegel_modeq_scalar**, **siegel_modeq_num** and
**siegel_modeq_den** fot Hilbert modular equations in Gundlach
invariants (both symmetric and nonsymmetric) in discriminant 5.

::
   slong hilbert_modeq_startprec(fmpq* params, slong ell, slong len)
   slong hilbert_modeq_nextprec(slong current_prec);

Manage complex precisions during the computations of Hilbert modular
equations.

::
   void hilbert_modeq_gundlach_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta);

Analogue of **siegel_modeq_rescale** in the case of Hilbert modular
equations in Gundlach invariants for discriminant 5.

::
   int hilbert_modeq_igusa_eval_Q(fmpz_poly_struct* num_vec,
			       fmpz_t den, fmpq* rs, slong ell, slong delta);
   int hilbert_modeq_nonsym_igusa_eval_Q(fmpz_poly_struct* num_vec,
				      fmpz_t den, fmpq* rs, slong ell, const fmpz_poly_t beta,
				      slong delta);
   int hilbert_modeq_gundlach_eval_Q(fmpz_poly_struct* num_vec,
				  fmpz_t den, fmpq* g, slong ell, slong delta);
   int hilbert_modeq_nonsym_gundlach_eval_Q(fmpz_poly_struct* num_vec,
					 fmpz_t den, fmpq* mn, slong ell,
					 const fmpz_poly_t beta, slong delta);
   int hilbert_modeq_igusa_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				const fmpz* rs, slong ell, slong delta,
				const fmpz_mod_ctx_t ctx);
   int hilbert_modeq_nonsym_igusa_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				       const fmpz* rs, slong ell, const fmpz_poly_t beta,
				       slong delta, const fmpz_mod_ctx_t ctx);
   int hilbert_modeq_gundlach_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				   const fmpz* g, slong ell,
				   slong delta, const fmpz_mod_ctx_t ctx);
   int hilbert_modeq_nonsym_gundlach_eval_Fp(fmpz_mod_poly_struct* pol_vec,
					  fmpq* mn, slong ell,
					  const fmpz_poly_t beta, slong delta,
					  const fmpz_mod_ctx_t ctx);

Evaluate Hilbert modular equations of various types.


