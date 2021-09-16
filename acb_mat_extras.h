
#ifndef ACB_MAT_EXTRAS_H
#define ACB_MAT_EXTRAS_H

#include "arf.h"
#include "arb.h"
#include "arb_mat.h"
#include "acb_mat.h"

void
acb_mat_get_real(arb_mat_t re, const acb_mat_t z);

void
acb_mat_get_imag(arb_mat_t im, const acb_mat_t z);

void
acb_mat_set_arb_arb(acb_mat_t z, const arb_mat_t re, const arb_mat_t im);

void
arb_mat_randtest_precise(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits);

void
arb_mat_randtest_sym_precise(arb_mat_t r, flint_rand_t state, slong prec, slong mag_bits);

void
arb_mat_randtest_sym_pos(arb_mat_t r, flint_rand_t state, slong prec);

int
arb_mat_is_nonsymmetric(const arb_mat_t m);

void
arb_mat_congr_fmpz_mat(arb_mat_t r, const fmpz_mat_t u, const arb_mat_t m, slong prec);

int
arb_mat_is_minkowski_reduced(const arb_mat_t r, const arb_t tol, slong prec);

int
arb_mat_not_minkowski_reduced(const arb_mat_t r, slong prec);

int
arb_mat_minkowski_reduce(arb_mat_t r, fmpz_mat_t u, const arb_mat_t m,
			 const arb_t tol, slong prec);

void
arb_mat_lambda(arb_t lambda, const arb_mat_t m, slong prec);

#endif

