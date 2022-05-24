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

#define MODEQ_CTX_ALLOC 100
#define MODEQ_MAX_NB_MONOMIALS 11
#define MODEQ_MAX_PREC n_pow(10,6)

#define SIEGEL_START_PREC_MUL 30
#define SIEGEL_START_PREC_ADD 0
#define SIEGEL_MUL_PREC 1.5

#define SIEGEL_2STEP_START_PREC_MUL 30
#define SIEGEL_2STEP_START_PREC_ADD 0

#define HILBERT_START_PREC_MUL 10
#define HILBERT_START_PREC_ADD 200
#define HILBERT_MUL_PREC 1.8


/* Define a data structure modeq_t that contains all the necessary
   data to translate a Hecke data structure into a bunch of polynomials */

typedef struct
{
  slong wt; /* Weight of Siegel modular forms used as projective coordinates */
  slong nb; /* Number of monomials */
  fmpz_mpoly_struct* monomials; /* Sufficient number of monomials of weight wt */
  fmpz_mpoly_t den; /* Denominator for the chosen absolute invariants */
  fmpz_mpoly_t num; /* Write all monomials in terms of parameter num/den */
  fmpz_mpoly_ctx_t ctx;
  slong nb_pairs; /* Number of pairs of maybe-equal projective coordinates */
  size_t alloc_pairs;
  slong* pairs;
} modeq_ctx_struct;

typedef modeq_ctx_struct modeq_ctx_t[1];

typedef struct
{
  slong d; /* Degree of Hecke correspondence */
  slong nb; /* Number of monomials */
  acb_poly_struct* all_nums; /* Equation satisfied by parameter,
			       and numerators of F where monomials = F(parameter) */
  acb_t den; /* Common denominator for all of these polynomials */
} modeq_acb_struct;

typedef modeq_acb_struct modeq_acb_t[1];

typedef struct
{
  slong d;
  slong nb;
  fmpz_poly_struct* all_nums;
  fmpz_t den;
} modeq_struct;

typedef modeq_struct modeq_t[1];


/* Handle modeq_ctx_t structures */

void modeq_ctx_init(modeq_ctx_t ctx);

void modeq_ctx_clear(modeq_ctx_t ctx);

#define modeq_ctx_weight(ctx) ((ctx)->wt)
#define modeq_ctx_nb(ctx) ((ctx)->nb)
#define modeq_ctx_monomial(ctx, k) (&(ctx)->monomials[(k)])
#define modeq_ctx_den(ctx) ((ctx)->den)
#define modeq_ctx_num(ctx) ((ctx)->num)
#define modeq_ctx_ctx(M) ((M)->ctx)
#define modeq_ctx_nb_pairs(ctx) ((ctx)->nb_pairs)
#define modeq_ctx_pair(ctx, k) (&(ctx)->pairs[2*(k)])
#define modeq_ctx_alloc_pairs(ctx) ((ctx)->alloc_pairs)

void modeq_ctx_add_pair(modeq_ctx_t ctx, slong i1, slong i2);

int modeq_ctx_is_pair(slong i1, slong i2, const modeq_ctx_t ctx);

int modeq_ctx_choose(modeq_ctx_t ctx, acb_srcptr I, slong nb, slong prec);

/* Handle modeq_t structures */

void modeq_init(modeq_t E);

void modeq_clear(modeq_t E);

void modeq_acb_init(modeq_acb_t E);

void modeq_acb_clear(modeq_acb_t E);

void modeq_acb_set(modeq_acb_t R, const modeq_acb_t E);

#define modeq_degree(E) ((E)->d)
#define modeq_nb(E) ((E)->nb)
#define modeq_all_nums(E) (E->all_nums)
#define modeq_equation(E) (&(E)->all_nums[0])
#define modeq_interpolate(E, j) (&(E)->all_nums[(j+1)])
#define modeq_den(E) ((E)->den)


/* Generic computations */

void modeq_product_trees(modeq_acb_t E, const hecke_t H,
			 const modeq_ctx_t ctx, slong prec);

void modeq_rescale(modeq_acb_t R, const modeq_acb_t E,
		   const acb_t c, slong prec);

void modeq_verbose_start(slong prec);

int modeq_stop(int res, slong prec);

int modeq_round(modeq_t R, const modeq_acb_t E);

int modeq_rationalize(modeq_t R, const modeq_acb_t E, slong prec);


void modeq_isog_monomials_Q(fmpz* M, const modeq_t E, const fmpq_t root,
			    slong mult);

void modeq_isog_monomials_Fp(fmpz* M, const modeq_t E, const fmpz_t root,
			     slong mult, const fmpz_mod_ctx_t ctx);

void modeq_isog_monomials_nf(fmpz_poly_struct* M, const modeq_t E,
			     const fmpz_poly_t factor, slong mult);

void modeq_all_isog_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
		      slong* exp_array, const modeq_t E, const modeq_ctx_t ctx);

int modeq_all_isog_Fp(slong* nb_roots, fmpz* all_M, slong* nb_M,
		      slong* exp_array, const modeq_t E, const modeq_ctx_t ctx,
		      const fmpz_mod_ctx_t fpctx);


/* Siegel modular equations */

slong siegel_modeq_startprec(fmpz* I, slong ell);

slong siegel_modeq_nextprec(slong current_prec);

void siegel_modeq_scalar(acb_t c, const hecke_t H, fmpz* I,
			 const modeq_ctx_t ctx, slong prec);

int siegel_modeq_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I, slong ell);

int siegel_isog_monomials_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
			    slong* exp_array, fmpz* I, slong ell);

int siegel_isog_monomials_Fp(slong* nb_roots, fmpz* all_M, slong* nb_M,
			     slong* exp_array, fmpz* I, slong ell,
			     const fmpz_mod_ctx_t ctx);


/* 2-step Siegel modular equations */

slong siegel_modeq_2step_startprec(fmpz* I, slong ell);

void siegel_modeq_2step_scalar(acb_t c, const hecke_t H, fmpz* I,
			       const modeq_ctx_t ctx, slong prec);

int siegel_modeq_2step_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I, slong ell);

int siegel_2step_isog_monomials_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
				  slong* exp_array, fmpz* I, slong ell);

int siegel_2step_isog_monomials_Fp(slong* nb_roots, fmpz* all_M, slong* nb_M,
                                   slong* exp_array, fmpz* I, slong ell,
				   const fmpz_mod_ctx_t ctx);


/* Hilbert modular equations */

slong hilbert_modeq_startprec(fmpz* I, slong ell, slong delta);

slong hilbert_modeq_nextprec(slong current_prec);

void hilbert_modeq_scalar(acb_t c, const hecke_t H, fmpz* I,
			  const modeq_ctx_t ctx, slong delta, slong prec);

int hilbert_modeq_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I,
		       slong ell, slong delta);

int hilbert_isog_monomials_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
			     slong* exp_array, fmpz* I, slong ell,
			     slong delta);

int hilbert_isog_monomials_Fp(slong* nb_roots, fmpz* all_M, slong* nb_M,
			      slong* exp_array, fmpz* I, slong ell,
			      slong delta, const fmpz_mod_ctx_t fpctx);


int hilbert_modeq_eval_split(modeq_t R1, modeq_t R2, modeq_ctx_t ctx1,
			     modeq_ctx_t ctx2, fmpz* I, slong ell, slong delta);


#endif 
