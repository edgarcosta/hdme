
#ifndef DERIV_H
#define DERIV_H

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
#include "modular.h"


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
