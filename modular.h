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

#include "hilbert.h"
#include "igusa.h"
#include "siegel.h"
#include "theta.h"

#define MODEQ_VERBOSE 1
#define MODEQ_RED_TOL_BITS 50
#define MODEQ_MAX_PREC n_pow(10,6)

#define SIEGEL_START_PREC_MUL 30
#define SIEGEL_START_PREC_ADD 0
#define SIEGEL_MUL_PREC 1.5

#define HILBERT_START_PREC_MUL 10
#define HILBERT_START_PREC_ADD 200
#define HILBERT_MUL_PREC 1.8


/* Generic functions for all types of modular equations */

void product_tree_1(acb_poly_t P, acb_srcptr xi, acb_srcptr yi, slong d, slong prec);

void product_tree_2(acb_poly_t Q, acb_srcptr xi, acb_srcptr yi, acb_srcptr zi,
		    slong d, slong prec);

int modeq_round_coeff(fmpz_t c, const acb_t x);

int modeq_round_poly(fmpz_poly_t pol, arf_t max_radius,
		     const acb_poly_t pol_acb, slong degree);

void modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong nb, slong prec);

int modeq_round(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* num_vec_acb,
		const acb_t den_acb, slong degree, slong nb);

int modeq_rational_coeff(fmpq_t c, fmpz_t den, const acb_t x,
			 const fmpz_t probable_den, slong prec);

int modeq_rational_poly(fmpq_poly_t pol, const acb_poly_t pol_acb,
			slong degree, slong prec);

int modeq_rational(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* pol_vec_acb,
		   slong degree, slong nb, slong prec);

void modeq_input_get_fmpz(fmpz_t den, fmpz* num, fmpq* j, slong len);

slong modeq_height_fmpz(const fmpz* j, slong len);

slong modeq_height_fmpq(fmpq* j, slong len);

void modeq_simplify(fmpz_poly_struct* num_vec, fmpz_t den, slong degree, slong nb);

void modeq_factor_Q(slong* nb_factors, fmpz_poly_struct* factors, slong* exps,
		    const fmpz_poly_t pol);

void modeq_roots_Q(slong* nb_roots, fmpq* roots, slong* mults,
		   const fmpz_poly_t pol);

int modeq_isog_invariants_Q(fmpq* j, const fmpz_poly_struct* num_vec,
			    const fmpq_t root, slong nb);

void modeq_factor_Fp(slong* nb_factors, fmpz_mod_poly_struct* factors, slong* exps,
		     const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

void modeq_roots_Fp(slong* nb_roots, fmpz* roots, slong* mults,
		    const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

int modeq_isog_invariants_Fp(fmpz* j, const fmpz_mod_poly_struct* pol_vec,
			     const fmpz_t root, slong nb,
			     const fmpz_mod_ctx_t ctx);

void modeq_input_lift(fmpq* j, const fmpz* input, slong nb);

int modeq_reduce(fmpz_mod_poly_struct* red_vec, const fmpz_poly_struct* num_vec,
		 const fmpz_t den, slong nb, const fmpz_mod_ctx_t ctx);

int modeq_isog_invariants_C(acb_ptr j, const fmpz_poly_struct* num_vec,
			    const acb_t root, slong nb, slong prec);

int modeq_isog_invariants_nf(fmpq_poly_struct* j, const fmpz_poly_struct* num_vec,
			     slong nb, const fmpz_poly_t field);

/* Siegel modular equations */

slong siegel_modeq_height_fmpz(const fmpz* j);

slong siegel_modeq_height_fmpq(fmpq* j);

slong siegel_nb_cosets(slong ell);

void siegel_coset(fmpz_mat_t m, slong k, slong ell);

int siegel_modeq_theta2(acb_ptr th2_vec, acb_ptr stardets,
			const acb_mat_t tau, slong ell, slong prec);

void siegel_modeq_exps(slong* e, slong* a, slong* b, slong ell);

void siegel_modeq_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
			 slong ell, slong prec);

void siegel_modeq_num(acb_poly_struct* num_vec_acb,
		      acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);

void siegel_modeq_den(acb_t den, acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);

slong siegel_modeq_startprec_fmpz(const fmpz* j, slong ell);

slong siegel_modeq_startprec_fmpq(fmpq* j, slong ell);

slong siegel_modeq_nextprec(slong current_prec);

void siegel_modeq_rescale(fmpz_t scal, fmpq* j, slong ell);

int siegel_modeq_eval_Q(fmpz_poly_struct* num_vec,
			fmpz_t den, fmpq* j, slong ell);

int siegel_modeq_eval_Fp(fmpz_mod_poly_struct* pol_vec,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx);

int siegel_modeq_eval_C(acb_poly_struct* num_vec, acb_t den, acb_srcptr j, slong ell,
			slong prec);			

int siegel_modeq_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j,
				   fmpq* j, slong ell);

int siegel_modeq_2step_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j,
					 fmpq* j, slong ell);

/* Hilbert modular equations */

slong hilbert_nb_cosets(slong ell, slong delta);

void hilbert_coset(fmpz_poly_mat_t m, slong k, slong ell, slong delta);

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



#endif 
