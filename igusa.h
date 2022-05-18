/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef IGUSA_H
#define IGUSA_H

#include <acb.h>
#include <acb_mat.h>
#include <acb_poly.h>
#include <arb.h>
#include <flint/fmpq_mpoly.h>
#include <flint/ulong_extras.h>

#include "hdme_data.h"
#include "siegel.h"
#include "theta.h"

#define COV_WEIGHTS {4,6,10,12}
#define THOMAE_LOWPREC 50
#define THOMAE_MULPREC 8
#define THOMAE_VERBOSE 0
/* #define COV_NB 4 would clutter the code. */


/* Igusa covariants from theta constants */

void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec);

void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec);

void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec);

void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec);

void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec);

#define cov_I4(I) &(I)[0]
#define cov_I6prime(I) &(I)[1]
#define cov_I10(I) &(I)[2]
#define cov_I12(I) &(I)[3]

void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec);

int cov_from_tau(acb_ptr I, const acb_mat_t tau, slong prec);


/* Weighted polynomials */

void cov_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx);

void cov_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx);

void cov_mpoly_print(const fmpz_mpoly_t pol, const fmpz_mpoly_ctx_t ctx);

void cov_monomial_exps(slong* exps, const fmpz_mpoly_t mon,
		       const fmpz_mpoly_ctx_t ctx);

void cov_monomial(fmpz_mpoly_t mon, slong* exps, const fmpz_mpoly_ctx_t ctx);

slong cov_nb_base_monomials(slong wt);

void cov_base_exps(slong* exps, slong wt, slong k);

void cov_base_monomial(fmpz_mpoly_t mon, slong wt, slong k,
		       const fmpz_mpoly_ctx_t ctx);

void cov_mpoly_eval(acb_t ev, const fmpz_mpoly_t pol, acb_srcptr I,
		    const fmpz_mpoly_ctx_t ctx, slong prec);

void cov_eval_base_monomials(acb_ptr ev, acb_srcptr I, slong wt, slong prec);


/* Rescaling */

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal, slong prec);

void cov_rescale_fmpz(fmpz* I, fmpz* S, const fmpz_t scal);

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal);

void cov_divexact_fmpz(fmpz* I, fmpz* S, const fmpz_t scal);

void cov_rescale_fmpz_si(fmpz* I, fmpz* S, slong scal);

void cov_divexact_fmpz_si(fmpz* I, fmpz* S, slong scal);

void cov_normalize_fmpz(fmpz* I, fmpz* S);

slong cov_height(const fmpz* I, slong len, slong* weights);

void cov_min_weight_combination(slong* wt, slong* i1, slong* i2,
				slong* e1, slong* e2, fmpz* I);

int cov_find_rescaling(acb_t scal, acb_srcptr I, fmpz* S, slong prec);

int cov_no_rescale_to_one(acb_srcptr I, slong prec);

int cov_distinct(acb_srcptr I1, acb_srcptr I2, slong prec);


/* Different covariants: classical Igusa--Clebsch I2, I4, I6, I10
   and Clebsch A, B, C, D */

void igusa_I2(acb_t I2, acb_srcptr I, slong prec);

int igusa_I2_fmpz(fmpz_t I2, fmpz* I);

void igusa_I6(acb_t I6, acb_srcptr I, slong prec);

int igusa_I6_fmpz(fmpz_t I6, fmpz* I);

void igusa_IC(acb_ptr IC, acb_srcptr I, slong prec);

void igusa_IC_fmpz(fmpz* IC, fmpz* I);

void igusa_I6prime(acb_t I6prime, acb_srcptr IC, slong prec);

int igusa_I6prime_fmpz(fmpz_t I6prime, fmpz* IC);

void cov_from_IC(acb_ptr I, acb_srcptr IC, slong prec);

void cov_from_IC_fmpz(fmpz* I, fmpz* IC);

void igusa_ABCD_from_IC(acb_ptr ABCD, acb_srcptr IC, slong prec);

void igusa_R2_from_IC(acb_t res, acb_srcptr IC, slong prec);

void igusa_ABCD_from_IC_fmpz(fmpq* ABCD, fmpz* IC);

void igusa_R2_from_IC_fmpz(fmpq_t R2, fmpz* IC);


/* Covariants from curve coefficients */

void curve_coeffs(acb_ptr ai, const acb_poly_t crv);

void curve_coeffs_fmpz(fmpz* ai, const fmpz_poly_t crv);

void cov_from_curve(acb_ptr I, const acb_poly_t crv, slong prec);

void cov_from_curve_fmpz(fmpz* I, const fmpz_poly_t crv);


/* Mestre's algorithm: use I2, I4, I6, I10 */

int igusa_has_generic_automorphisms(acb_srcptr IC, slong prec);

void igusa_generic_randtest(acb_poly_t crv, acb_ptr IC, flint_rand_t state, slong prec);

void mestre_conic(acb_ptr conic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

void mestre_conic_randtest(acb_ptr conic, flint_rand_t state, slong prec);

void mestre_cubic(acb_ptr cubic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec);

void mestre_subst_in_conic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr conic, slong prec);

void mestre_subst_in_cubic(acb_poly_t subst, const acb_poly_t x1, const acb_poly_t x2,
			   const acb_poly_t x3, acb_srcptr cubic, slong prec);

int mestre_point_on_conic(acb_ptr pt, acb_srcptr conic, slong prec);

int mestre_point_is_outside_conic(acb_srcptr pt, acb_srcptr conic, slong prec);

void mestre_eval_cubic(acb_t res, acb_srcptr pt, acb_srcptr cubic, slong prec);

void mestre_parametrize_conic(acb_poly_t x1, acb_poly_t x2, acb_poly_t x3,
			     acb_srcptr pt, acb_srcptr conic, slong prec);

int mestre(acb_poly_t crv, acb_srcptr IC, slong prec);

void cardona(acb_poly_t crv, acb_srcptr IC, slong prec);


/* Thomae's formulae: back to I4, I6prime, I10, I12 */

int thomae_roots(acb_ptr roots, const acb_poly_t crv, slong prec);

void thomae_reorder(acb_ptr new_roots, acb_srcptr roots, slong perm);

void thomae_rosenhain(acb_ptr ros, acb_srcptr roots, slong prec);

void thomae_theta4(acb_ptr th4, acb_srcptr ros, slong prec);

void thomae_theta2(acb_ptr th2, acb_srcptr th4, acb_srcptr ros, slong signs, slong prec);

int thomae_discard(acb_srcptr th2, slong prec);

int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec);

slong thomae_startprec(slong prec);

int thomae_correct_signs(slong* perm, slong* signs, acb_srcptr roots,
			 acb_srcptr I, slong prec);

int tau_theta2_from_curve(acb_mat_t tau, acb_ptr th2, const acb_poly_t crv,
			  slong prec);

int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec);

int tau_theta2_from_igusa(acb_mat_t tau, acb_ptr th2, acb_srcptr I, slong prec);


/* Period computations for products of elliptic curves, and more
   generally abelian surfaces with extra automorphisms */

int cov_is_g2_curve(acb_srcptr I);

int cov_is_g2_curve_fmpz(fmpz* I);

void igusa_ec_j1j2(acb_ptr j, fmpz* I, slong prec);

int igusa_ec_period(acb_t tau, const acb_t j, slong prec);

int tau_theta2_from_igusa_ec(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec);

int mestre_fmpz(acb_poly_t crv, fmpz* IC, slong prec);

int tau_theta2_from_igusa_fmpz(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec);


#endif 
