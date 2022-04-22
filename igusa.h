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

#define THOMAE_LOWPREC 50
#define THOMAE_MULPREC 8
#define THOMAE_VERBOSE 0


/* Igusa covariants from theta constants */

void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec);

void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec);

void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec);

void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec);

void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec);

void igusa_h(acb_ptr h, acb_srcptr theta2, slong prec);

void cov_from_h(acb_ptr I, acb_srcptr h, slong prec);

void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec);

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal, slong prec);

void igusa_from_theta2(acb_ptr j, acb_srcptr theta2, slong prec);

int igusa_from_tau(acb_ptr j, const acb_mat_t tau, slong prec);

int igusa_is_defined(acb_srcptr j);


/* Igusa covariants as integers */

void igusa_from_cov(acb_ptr j, acb_srcptr I, slong prec);

void igusa_from_cov_fmpz(fmpq* j, const fmpz* I);

void cov_from_igusa(acb_ptr I, acb_srcptr j, slong prec);

void cov_rescale_fmpz(fmpz* I, fmpz* S, const fmpz_t scal);

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal);

void cov_divexact_fmpz(fmpz* I, fmpz* S, const fmpz_t scal);

void cov_rescale_fmpz_si(fmpz* I, fmpz* S, slong scal);

void cov_divexact_fmpz_si(fmpz* I, fmpz* S, slong scal);

void cov_normalize_fmpz(fmpz* I, fmpz* S);

void cov_from_igusa_fmpz(fmpz* I, fmpq* j);


/* Igusa covariants from curve coefficients */

void curve_coeffs(acb_ptr ai, const acb_poly_t crv);

void curve_coeffs_fmpz(fmpz* ai, const fmpz_poly_t crv);

void igusa_scalar_covariants(acb_ptr I, const acb_poly_t crv, slong prec);

void igusa_scalar_covariants_fmpz(fmpz* I, const fmpz_poly_t crv);

void igusa_from_curve(acb_ptr j, const acb_poly_t crv, slong prec);

void igusa_from_curve_fmpz(fmpq* j, const fmpz_poly_t crv);


/* Different covariants: I6, and Clebsch */

void igusa_I6(acb_t I6, acb_srcptr I, slong prec);

int igusa_I6_fmpz(fmpz_t I6, fmpz* I);

void igusa_switch_I6_fmpz(fmpz* I, fmpz* S);

void igusa_I6prime(acb_t I6prime, acb_srcptr I, slong prec);

int igusa_I6prime_fmpz(fmpz_t I6prime, fmpz* I);

void igusa_switch_I6prime_fmpz(fmpz* I, fmpz* S);

void igusa_clebsch(acb_ptr ABCD, acb_srcptr I, slong prec);

void igusa_R2(acb_t res, acb_srcptr I, slong prec);


/* Mestre's algorithm */

int igusa_has_generic_automorphisms(acb_srcptr I, slong prec);

void igusa_generic_randtest(acb_poly_t crv, acb_ptr I, flint_rand_t state, slong prec);

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

int mestre(acb_poly_t crv, acb_srcptr I, slong prec);


/* Thomae's formulae */

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

int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec);

int tau_theta2_from_igusa(acb_mat_t tau, acb_ptr th2, acb_srcptr I, slong prec);


#endif 
