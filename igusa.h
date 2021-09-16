
#ifndef IGUSA_H
#define IGUSA_H

#include "ulong_extras.h"
#include "arb.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_mat_extras.h"
#include "siegel.h"
#include "theta.h"

#define THOMAE_LOWPREC 50
#define THOMAE_MULPREC 4
#define THOMAE_VERBOSE 0

/* Igusa invariants from theta constants */

void igusa_h4(acb_t h4, acb_srcptr theta2, slong prec);

void igusa_h6(acb_t h6, acb_srcptr theta2, slong prec);

void igusa_h10(acb_t h10, acb_srcptr theta2, slong prec);

void igusa_h12(acb_t h12, acb_srcptr theta2, slong prec);

void igusa_h16(acb_t h16, acb_srcptr theta2, slong prec);

void igusa_h(acb_ptr h, acb_srcptr theta2, slong prec);

void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec);

void igusa_from_theta2(acb_ptr j, acb_srcptr theta2, slong prec);

int igusa_from_tau(acb_ptr j, const acb_mat_t tau, slong prec);

int igusa_is_defined(acb_srcptr j);

/* Igusa covariants from curve coefficients */

void igusa_I2_autogen(acb_t res, const acb_t a0, const acb_t a1,
		      const  acb_t a2, const acb_t a3, const acb_t a4,
		      const acb_t a5, const acb_t a6, slong prec);

void igusa_I4_autogen(acb_t res, const acb_t a0, const acb_t a1,
		      const  acb_t a2, const acb_t a3, const acb_t a4,
		      const acb_t a5, const acb_t a6, slong prec);

void igusa_I6prime_autogen(acb_t res, const acb_t a0, const acb_t a1,
			   const  acb_t a2, const acb_t a3, const acb_t a4,
			   const acb_t a5, const acb_t a6, slong prec);

void igusa_I10_autogen(acb_t res, const acb_t a0, const acb_t a1,
		       const  acb_t a2, const acb_t a3, const acb_t a4,
		       const acb_t a5, const acb_t a6, slong prec);

void igusa_I2_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1,
		   const  fmpz_t a2, const fmpz_t a3, const fmpz_t a4,
		   const fmpz_t a5, const fmpz_t a6);

void igusa_I4_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1,
		   const  fmpz_t a2, const fmpz_t a3, const fmpz_t a4,
		   const fmpz_t a5, const fmpz_t a6);

void igusa_I6prime_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1,
			const  fmpz_t a2, const fmpz_t a3, const fmpz_t a4,
			const fmpz_t a5, const fmpz_t a6);

void igusa_I10_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1,
		    const  fmpz_t a2, const fmpz_t a3, const fmpz_t a4,
		    const fmpz_t a5, const fmpz_t a6);

void igusa_scalar_covariants(acb_ptr I, const acb_poly_t crv, slong prec);

void igusa_scalar_covariants_fmpz(fmpz* I, const fmpz_poly_t crv);

void igusa_from_cov(acb_ptr j, acb_srcptr I, slong prec);

void igusa_from_cov_fmpz(fmpq* j, const fmpz* I);

void cov_from_igusa(acb_ptr I, acb_srcptr j, slong prec);

void igusa_from_curve(acb_ptr j, const acb_poly_t crv, slong prec);

void igusa_from_curve_fmpz(fmpq* j, const fmpz_poly_t crv);

/* Different covariants: I6, and Clebsch */

void igusa_I6(acb_t I6, acb_srcptr I, slong prec);

void igusa_clebsch(acb_ptr ABCD, acb_srcptr I, slong prec);

void igusa_R2_autogen(acb_t res, const acb_t I2, const acb_t I4, const acb_t I6,
		      const acb_t I10, slong prec);

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
