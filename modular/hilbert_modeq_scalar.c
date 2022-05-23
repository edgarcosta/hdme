
#include "modular.h"

/* See also siegel_modeq_scalar */

void hilbert_modeq_scalar(acb_t c, const hecke_t H, fmpz* I,
			  const modeq_ctx_t ctx, slong delta, slong prec)
{
  slong wt = modeq_ctx_weight(ctx);
  slong ell = hecke_ell(H);
  fmpz* G;
  acb_ptr G_tau;
  slong weights[3] = GUNDLACH_WEIGHTS_5;
  acb_t s;
  slong k;
  int res;

  G = _fmpz_vec_init(3);
  G_tau = _acb_vec_init(3);
  acb_init(res);

  if (delta != 5)
    {
      flint_printf("(hilbert_modeq_scalar) Error: delta=5 is required\n");
      fflush(stdout);
      flint_abort();
    }
  
  gundlach_from_igusa(G_tau, hecke_I_tau(H), delta, prec);
  gundlach_from_igusa_fmpz(G, I, delta, prec);
  
  /* Step 1: product of all stardets to the correct weight, to get
     a Hilbert modular form */
  acb_one(c);
  for (k = 0; k < hecke_nb(H); k++)
    {
      acb_inv(s, hecke_stardet(H, k), prec);
      acb_mul_si(s, s, ell, prec);
      acb_mul(c, c, s, prec);
    }
  
  /* Step 2: adjust power of ell that keeps integer coefficients */
  acb_set_si(s, ell);
  acb_pow_si(s, s, -2);
  acb_mul(c, c, s, prec);

  
  /* Step 3: multiply by scaling factor between hecke_I_tau and I */
  cov_find_rescaling(s, G_tau, G, 3, weights, prec);  
  acb_mul(c, c, s, prec);
  
  acb_pow_si(c, c, wt, prec);

  _fmpz_vec_clear(G, 3);
  _acb_vec_clear(G_tau, 3);
  acb_clear(s);
}
