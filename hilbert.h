
#ifndef HILBERT_H
#define HILBERT_H

#include <stdio.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz_poly_mat.h"
#include "fmpq_mpoly.h"
#include "fmpz_lll.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_mat_extras.h"
#include "siegel.h"
#include "theta.h"
#include "igusa.h"

#define HILBERT_MAX_STRLEN 4096
#define HUMBERT_DATA_PATH HDME_PATH"/data/humbert"
#define HILBERT_DATA_PATH HDME_PATH"/data/hilbert"
#define GUNDLACH_DATA_PATH HDME_PATH"/data/gundlach"

int hilbert_is_fundamental(slong delta);

int hilbert_is_totally_positive(const fmpz_poly_t x, slong delta);

int hilbert_splits(fmpz_poly_t beta, slong ell, slong delta);

void hilbert_conjugate(fmpz_poly_t xbar, fmpz_poly_t x, slong delta);

char** humbert_vars_init();

void humbert_vars_clear(char** vars);

void humbert_vars_set(char** vars, slong delta);

void humbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
		       slong delta, const fmpq_mpoly_ctx_t ctx);

void hilbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
		       slong delta, const fmpq_mpoly_ctx_t ctx);

void gundlach_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
			slong delta, const fmpq_mpoly_ctx_t ctx);

void fmpq_mpoly_evaluate_all_acb(acb_t ev, const fmpq_mpoly_t pol, acb_srcptr vals,
				 const fmpq_mpoly_ctx_t ctx, slong prec);

void humbert_AA1BB1B2(acb_ptr AA1BB1B2, const acb_t r, const acb_t s, slong delta,
		      slong prec);

void humbert_cov_from_AA1BB1B2(acb_ptr I, acb_srcptr AA1BB1B2, slong prec);

void humbert_parametrize(acb_ptr I, const acb_t r, const acb_t s, slong delta,
			 slong prec);

void hilbert_parametrize(acb_ptr I, const acb_t r, const acb_t s, slong delta, slong prec);

void gundlach_from_igusa(acb_ptr g, acb_srcptr I, slong delta, slong prec);

void igusa_from_gundlach(acb_ptr j, acb_srcptr g, slong delta, slong prec);

void hilbert_halfspace_randtest(acb_t t1, acb_t t2, flint_rand_t state, slong prec);

void hilbert_R(acb_mat_t R, slong delta, slong prec);

void hilbert_map(acb_mat_t tau, const acb_t t1, const acb_t t2, slong delta, slong prec);

int hilbert_linear_combination(fmpz* abcde, const acb_mat_t tau, slong delta, slong prec);

int hilbert_inverse(acb_t t1, acb_t t2, sp2gz_t eta, const acb_mat_t tau,
		    slong delta, slong prec);

void hilbert_sigma1(acb_t z, const fmpz_poly_t x, slong delta, slong prec);

void hilbert_sigma2(acb_t z, const fmpz_poly_t x, slong delta, slong prec);

void hilbert_star(acb_t z, const fmpz_poly_mat_t m, const acb_t t1, const acb_t t2,
		  slong delta, slong prec);

void hilbert_transform(acb_t z1, acb_t z2, const fmpz_poly_mat_t m, const acb_t t1,
		       const acb_t t2, slong delta, slong prec);

void hilbert_scalar_mul(acb_t z1, acb_t z2, const fmpz_poly_t x, const acb_t t1,
			const acb_t t2, slong delta, slong prec);

void hilbert_scalar_div(acb_t z1, acb_t z2, const fmpz_poly_t x, const acb_t t1,
			const acb_t t2, slong delta, slong prec);

#endif
