
#ifndef MODULAR_H
#define MODULAR_H

#include "flint.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_poly_mat.h"
#include "ulong_extras.h"
#include "arb.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_mat_extras.h"
#include "siegel.h"
#include "theta.h"
#include "igusa.h"
#include "hilbert.h"

#define SIEGEL_START_PREC_MUL 30
#define SIEGEL_START_PREC_ADD 0
#define SIEGEL_MUL_PREC 1.5
#define SIEGEL_RED_TOL_BITS 50
#define SIEGEL_MAX_PREC n_pow(10,6)
/* Todo: use this in functions. */
#define SIEGEL_VERBOSE 1 

#define HILBERT_START_PREC_MUL 20
#define HILBERT_START_PREC_ADD 150
#define HILBERT_MUL_PREC 1.8
#define HILBERT_MAX_PREC n_pow(10,6)
#define HILBERT_VERBOSE 1

/* Siegel modular equations */

slong siegel_nb_cosets(slong ell);

void siegel_coset(sp2gz_t m, slong k, slong ell);

int siegel_modeq_theta2(acb_ptr th2_vec, acb_ptr stardets,
			const acb_mat_t tau, slong ell, slong prec);

void siegel_modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong ell, slong prec);

void product_tree_1(acb_poly_t P, acb_srcptr xi, acb_srcptr yi, slong d, slong prec);

void product_tree_2(acb_poly_t Q, acb_srcptr xi, acb_srcptr yi, acb_srcptr zi,
		    slong d, slong prec);

void siegel_modeq_exps(slong* e, slong* a, slong* b, slong ell);

void siegel_modeq_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
			 slong ell, slong prec);

void siegel_modeq_num(acb_poly_t num1, acb_poly_t num2, acb_poly_t num3,
		      acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);

void siegel_modeq_den(acb_t den, acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec);

int siegel_modeq_round_coeff(fmpz_t c, const acb_t x);

int siegel_modeq_round_poly(fmpz_poly_t pol, arf_t max_radius,
			    const acb_poly_t pol_acb, slong degree);

int siegel_modeq_round(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
		       fmpz_t den, const acb_poly_t num1_acb, const acb_poly_t num2_acb,
		       const acb_poly_t num3_acb, const acb_t den_acb, slong ell);

void siegel_modeq_fmpz_input(fmpz_t den, fmpz* num, fmpq* j, slong len);

slong siegel_modeq_height_fmpz(const fmpz* j);

slong siegel_modeq_startprec_fmpz(const fmpz* j, slong ell);

slong siegel_modeq_height_fmpq(fmpq* j);

slong siegel_modeq_startprec_fmpq(fmpq* j, slong ell);

slong siegel_modeq_nextprec(slong current_prec);

void siegel_modeq_fmpq_rescale(fmpz_t scal, fmpq* j, slong ell);

void siegel_modeq_simplify(fmpz_poly_t num1, fmpz_poly_t num2,
			   fmpz_poly_t num3, fmpz_t den, slong ell);

int siegel_modeq_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
			 fmpz_t den, fmpq* j, slong ell);

void siegel_modeq_factor_Q(slong* nb_factors, fmpz_poly_struct* factors, slong* exps,
			   const fmpz_poly_t pol);

void siegel_modeq_roots_Q(slong* nb_roots, fmpq* roots, slong* mults,
			  const fmpz_poly_t pol);

int siegel_modeq_isog_igusa_Q(fmpq* j, const fmpz_poly_t num1, const fmpz_poly_t num2,
			      const fmpz_poly_t num3, const fmpq_t root);

int siegel_modeq_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2, fmpz_mod_poly_t pol3,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx);

void siegel_modeq_factor_Fp(slong* nb_factors, fmpz_mod_poly_struct* factors, slong* exps,
			    const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

void siegel_modeq_roots_Fp(slong* nb_roots, fmpz* roots, slong* mults,
			   const fmpz_mod_poly_t pol, const fmpz_mod_ctx_t ctx);

int siegel_modeq_isog_igusa_Fp(fmpz* j, const fmpz_mod_poly_t pol1, const fmpz_mod_poly_t pol2,
			       const fmpz_mod_poly_t pol3, const fmpz_t root,
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



/* Hilbert modular equations */

slong hilbert_nb_cosets(slong ell, slong delta);

void hilbert_coset(fmpz_poly_mat_t m, slong k, slong ell, slong delta);

int hilbert_modeq_theta2(acb_ptr th2_vec, const acb_t t1, const acb_t t2,
			 const fmpz_poly_t beta, slong ell, slong delta, slong prec);

int hilbert_modeq_theta2_star(acb_ptr th2_vec, acb_ptr stardets,
			      const acb_t t1, const acb_t t2,
			      const fmpz_poly_t beta, slong ell, slong delta, slong prec);

void hilbert_modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong ell,
		       slong delta, slong prec);

void hilbert_modeq_sym_igusa_C(acb_poly_t pol1, acb_poly_t pol2, acb_poly_t pol3,
			       acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
			       slong ell, slong delta, slong prec);

void hilbert_modeq_nonsym_igusa_C(acb_poly_t pol1, acb_poly_t pol2,
				  acb_poly_t pol3, acb_srcptr I_vec, slong ell,
				  slong delta, slong prec);

void hilbert_modeq_gundlach_exps(slong* e, slong* a, slong* b, slong ell, slong delta);

void hilbert_modeq_gundlach_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
				   slong ell, slong delta, slong prec);

void hilbert_modeq_gundlach_num(acb_poly_t num1, acb_poly_t num2,
				acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
				const acb_t scal,
				slong ell, slong delta, slong prec);

void hilbert_modeq_gundlach_den(acb_t den, acb_srcptr I_vec_beta,
				acb_srcptr I_vec_betabar, const acb_t scal,
				slong ell, slong delta, slong prec);

int hilbert_modeq_coeff_Q(fmpq_t c, fmpz_t den, const acb_t x,
			  const fmpz_t probable_den, slong prec);

int hilbert_modeq_poly_Q(fmpq_poly_t num, const acb_poly_t pol_acb,
			 slong degree, slong prec);

int hilbert_modeq_gundlach_round(fmpz_poly_t num1, fmpz_poly_t num2,
				 fmpz_t den, const acb_poly_t num1_acb,
				 const acb_poly_t num2_acb,
				 const acb_t den_acb, slong ell, slong delta);

slong hilbert_modeq_height(fmpq* params, slong len);

slong hilbert_modeq_sym_igusa_startprec(fmpq* params, slong ell, slong len);

slong hilbert_modeq_nextprec(slong current_prec);

void hilbert_modeq_gundlach_fmpq_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta);

void hilbert_modeq_gundlach_simplify(fmpz_poly_t num1, fmpz_poly_t num2,
				     fmpz_t den, slong ell, slong delta);

int hilbert_modeq_sym_igusa_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
				   fmpz_t den, fmpq* rs, slong ell, slong delta);

int hilbert_modeq_nonsym_igusa_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
				      fmpz_t den, fmpq* rs, slong ell, const fmpz_poly_t beta,
				      slong delta);

int hilbert_modeq_gundlach_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2,
				  fmpz_t den, fmpq* g, slong ell, slong delta);

int hilbert_modeq_sym_igusa_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2,
				    fmpz_mod_poly_t pol3,
				    const fmpz* rs, slong ell, slong delta,
				    const fmpz_mod_ctx_t ctx);

int hilbert_modeq_nonsym_igusa_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2,
				       fmpz_mod_poly_t pol3,
				       const fmpz* rs, slong ell, const fmpz_poly_t beta,
				       slong delta, const fmpz_mod_ctx_t ctx);

int hilbert_modeq_gundlach_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2,
				   const fmpz* g, slong ell,
				   slong delta, const fmpz_mod_ctx_t ctx);

/* Derivatives of Hilbert modular equations */
		   
#endif 
