/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef MODULAR_H
#define MODULAR_H

#include <acb.h>
#include <acb_mat.h>
#include <acb_poly.h>
#include <arb.h>
#include <flint/flint.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/ulong_extras.h>

#include "siegel.h"
#include "theta.h"
#include "igusa.h"
#include "hilbert.h"
#include "hecke.h"

/* one may want to define verbosity at compile time */
#ifndef MODEQ_VERBOSE
#define MODEQ_VERBOSE 0
#endif

#define MODEQ_MAX_NB_MONOMIALS 11

#define MODEQ_RED_TOL_BITS 50
#define MODEQ_MAX_PREC n_pow(10,6)

#define SIEGEL_START_PREC_MUL 30
#define SIEGEL_START_PREC_ADD 0
#define SIEGEL_MUL_PREC 1.5

#define HILBERT_START_PREC_MUL 10
#define HILBERT_START_PREC_ADD 200
#define HILBERT_MUL_PREC 1.8


/* Define a data structure modeq_t that contains all the necessary
   data to translate a Hecke data structure into a bunch of polynomials */

typedef struct
{
  slong wt; /* Weight of Siegel modular forms used as projective coordinates */
  fmpq_mpoly_ctx_t mpolyctx; /* Context for polynomials in Igusa-Clebsch covariants */
  slong nb_monomials;
  fmpq_mpoly_struct* monomials; /* Sufficient number of monomials of weight wt */
  fmpq_mpoly_t den; /* Denominator for the chosen absolute invariants */
  fmpq_mpoly_t num; /* Write all monomials in terms of parameter num/den */
} modeq_ctx_struct;

typedef modeq_ctx_struct modeq_ctx_t[1];

typedef struct
{  
  /* Complex-valued polynomials */
  slong d;
  acb_ptr param; /* Complex values of num/den at isogenous period matrices */
  acb_poly_t eq_acb; /* Equation satisfied by parameter */
  acb_poly_struct* interp_acb; /* Numerators of F where monomials = F(parameter) */
  acb_t den_acb; /* Common denominator for all of these polynomials */
  
  acb_t rescale; /* Scalar cofactor to make modeqs_acb & den_acb integral, if known */
  acb_poly_t eq_resc;
  acb_poly_struct* interp_resc;
  acb_t den_resc;

  /* Rational-valued polynomials, when appropriate */
  fmpz_poly_t eq;
  fmpz_poly_struct* interp;
  fmpz_t den;
} modeq_struct;

typedef modeq_struct modeq_t[1];


void modeq_input_lift(fmpq* j, const fmpz* input, slong nb);

/* Handle modeq_ctx_t structures */

void modeq_ctx_init(modeq_ctx_t ctx);

void modeq_ctx_clear(modeq_ctx_t ctx);

#define modeq_ctx_weight(ctx) ((ctx)->wt)
#define modeq_ctx_mpolyctx(ctx) ((ctx)->mpolyctx)
#define modeq_ctx_nb_monomials(ctx) ((ctx)->nb_monomials)
#define modeq_ctx_monomial(ctx, k) (&(ctx)->monomials[(k)])
#define modeq_ctx_den(ctx) ((ctx)->den)
#define modeq_ctx_num(ctx) ((ctx)->num)

int modeq_ctx_set(modeq_ctx_t ctx, acb_srcptr I, slong nb);


/* Handle modeq_t structures */

void modeq_init(modeq_t E);

void modeq_clear(modeq_t E);

#define modeq_degree(E) ((E)->d)
#define modeq_param(E, k) (&(E)->param[(k)])
#define modeq_eq_acb(E) ((E)->eq_acb)
#define modeq_interp_acb(E, j) (&(E)->interp_acb[(j)])
#define modeq_den_acb(E) ((E)->den_acb)
#define modeq_rescale(E) ((E)->rescale)
#define modeq_eq_resc(E) ((E)->eq_resc)
#define modeq_interp_resc(E, j) (&(E)->interp_resc[(j)])
#define modeq_den_resc(E) ((E)->den_resc)
#define modeq_equation(E) ((E)->eq)
#define modeq_interpolate(E, j) (&(E)->interp[(j)])
#define modeq_den(E) ((E)->den)


/* Generic computations */

void modeq_product_trees(modeq_t E, const hecke_t H,
			 const modeq_ctx_t ctx, slong prec);

int modeq_round(modeq_t E);

int modeq_rational(modeq_t E);

/* void modeq_rescale_fmpq(fmpz_t den, fmpz* num, fmpq* j, slong len); */

int modeq_isog_invariants_Q(fmpz* j, const modeq_t E, const fmpq_t root);

int modeq_isog_invariants_Fp(fmpz* j, const modeq_t E,
			     const fmpz_t root, const fmpz_mod_ctx_t ctx);

int modeq_isog_invariants_C(acb_ptr j, const modeq_t E,
			    const acb_t root, slong prec);

int modeq_isog_invariants_nf(fmpz_poly_struct* j, const modeq_t E,
			     const fmpz_poly_t field);


/* Siegel modular equations specifics */

slong siegel_modeq_startprec(const fmpz* I, slong ell);

slong siegel_modeq_nextprec(slong current_prec);

void siegel_modeq_rescale(modeq_t E, const hecke_t H, const modeq_ctx_t ctx,
			      slong prec);

void siegel_modeq_exps(slong* e, slong* a, slong* b, slong ell);

void siegel_modeq_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
			 slong ell, slong prec);


void siegel_modeq_rescale(fmpz_t scal, fmpq* j, slong ell);


int siegel_modeq_eval_Q(modeq_t E, fmpz* I, slong ell);

int siegel_modeq_eval_C(modeq_t E, acb_srcptr j, slong ell, slong prec);


int siegel_modeq_isog_invariants_Q(slong* nb_roots, slong* weight,
				   fmpz* all_isog_j, fmpz* I, slong ell);

int siegel_modeq_2step_isog_invariants_Q(slong* nb_roots, slong* weight,
					 fmpz* all_isog_j, fmpz* I, slong ell);

/* Hilbert modular equations */

int hilbert_modeq_theta2(acb_ptr th2_vec, acb_srcptr t,
			 const fmpz_poly_t beta, slong ell, slong delta, slong prec);

int hilbert_modeq_theta2_star(acb_ptr th2_vec, acb_ptr stardets,
			      acb_srcptr t,
			      const fmpz_poly_t beta, slong ell, slong delta, slong prec);

void hilbert_modeq_igusa_C(acb_poly_struct* pol_vec,
			   acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
			   slong ell, slong delta, slong prec);

void hilbert_modeq_nonsym_igusa_C(acb_poly_struct* pol_vec, acb_srcptr I_vec, slong ell,
				  slong delta, slong prec);

void hilbert_modeq_gundlach_exps(slong* e, slong* a, slong* b, slong ell, slong delta);

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

slong hilbert_modeq_startprec(fmpq* params, slong ell, slong len);

slong hilbert_modeq_nextprec(slong current_prec);

void hilbert_modeq_gundlach_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta);

void hilbert_modeq_nonsym_gundlach_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta);

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


/* Derivatives of Siegel modular equations */

int siegel_modeq_dtheta(acb_ptr thvec, acb_ptr thder, const acb_mat_t tau,
				   slong ell, slong prec);

int siegel_modeq_dcov(acb_ptr Ider, acb_srcptr thvec, acb_srcptr thder,
			 slong k, slong prec);

int siegel_modeq_dnum(acb_poly_struct* nums, const acb_t scal, acb_srcptr thvec,
		      acb_srcptr thder, slong ell, slong prec);

int siegel_modeq_dden(acb_t den, const acb_t scal, acb_srcptr thvec,
		      acb_srcptr thder, slong ell, slong prec);

int siegel_modeq_dround(fmpz_poly_struct* nums, fmpz_t den, const acb_poly_struct nums_acb,
			const acb_t den_acb, slong prec);

int siegel_modeq_deval_zz(fmpz_poly_struct* nums, fmpz_t den, const fmpq* j,
			   slong ell);

int siegel_modeq_deval_fp(fmpz_mod_poly_struct* pols, const fmpz* j, slong ell,
			  const fmpz_mod_ctx_t ctx);

/* Todo: also return evaluated modular equations? Do both versions? */



/* Derivatives of Hilbert modular equations */

/* OLD */

/* int modeq_round(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* num_vec_acb, */
/* 		const acb_t den_acb, slong degree, slong nb); */


/* void modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong nb, slong prec); */

/* int modeq_rational(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* pol_vec_acb, */
/* 		   slong degree, slong nb, slong prec); */


/* void modeq_input_get_fmpz(fmpz_t den, fmpz* num, fmpq* j, slong len); */

/* slong modeq_height_fmpq(fmpq* j, slong len); */




/* void modeq_input_lift(fmpq* j, const fmpz* input, slong nb); */

/* slong siegel_modeq_height_fmpz(const fmpz* j); */

/* slong siegel_modeq_height_fmpq(fmpq* j); */

/* int siegel_modeq_theta2(acb_ptr th2_vec, acb_ptr stardets, */
/* 			const acb_mat_t tau, slong ell, slong prec); */

/* void siegel_modeq_num(acb_poly_struct* num_vec_acb, */
/* 		      acb_srcptr I_vec, const acb_t scal, */
/* 		      slong ell, slong prec); */

/* void siegel_modeq_den(acb_t den, acb_srcptr I_vec, const acb_t scal, */
/* 		      slong ell, slong prec); */

/* slong siegel_modeq_startprec_fmpz(const fmpz* j, slong ell); */

/* slong siegel_modeq_startprec_fmpq(fmpq* j, slong ell); */

/* int siegel_modeq_eval_Q(fmpz_poly_struct* num_vec, */
/* 			fmpz_t den, fmpq* j, slong ell); */

/* int siegel_modeq_eval_Fp(fmpz_mod_poly_struct* pol_vec, */
/* 			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx); */

/* int siegel_modeq_eval_C(acb_poly_struct* num_vec, acb_t den, acb_srcptr j, slong ell, */
/* 			slong prec); */

/* int siegel_modeq_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j, */
/* 				   fmpq* j, slong ell); */

/* int siegel_modeq_2step_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j, */
/* 					 fmpq* j, slong ell); */

#endif 
