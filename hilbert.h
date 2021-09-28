
#ifndef HILBERT_H
#define HILBERT_H

#include <stdio.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpq_mpoly.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_mat_extras.h"
#include "siegel.h"
#include "theta.h"
#include "igusa.h"

#define HILBERT_MAX_STRLEN 4096
#define HILBERT_DATA_PATH HDME_PATH"/data/humbert"

int hilbert_is_fundamental(slong delta);

char** humbert_vars_init();

void humbert_vars_clear(char** vars);

void humbert_vars_set(char** vars, slong delta);

void humbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
		       slong delta, const fmpq_mpoly_ctx_t ctx);

void humbert_param(acb_ptr cov, const acb_t r, const acb_t s, slong Delta,
		   slong prec);

#endif
