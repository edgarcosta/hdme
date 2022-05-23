
#include "modular.h"

int modeq_all_isog_Fp(slong* nb_roots, fmpz* all_M, slong* nb_M,
		       slong* exp_array, const modeq_t E, const modeq_ctx_t ctx,
		       const fmpz_mod_ctx_t fpctx)
{
  slong nb = modeq_degree(E);
  slong* mults;
  fmpz* roots;
  fmpz_mod_poly_t red;
  slong k;
  int res;

  mults = flint_malloc(nb * sizeof(slong));
  roots = _fmpz_vec_init(nb);
  fmpz_mod_poly_init(red, fpctx);
  
  res = pol_reduce_Fp(red, modeq_equation(E), modeq_den(E), fpctx);
  if (res)
    {
      pol_roots_Fp(nb_roots, roots, mults, red, fpctx);
      *nb_M = modeq_ctx_nb(ctx);
      for (k = 0; k < *nb_M; k++)
	{
	  cov_monomial_degrees(&exp_array[4*k], modeq_ctx_monomial(ctx, k),
			       modeq_ctx_ctx(ctx));
	}
      for (k = 0; k < *nb_roots; k++)
	{
	  modeq_isog_monomials_Fp(&all_M[(*nb_M) * k], E, &roots[k], mults[k]);
	}
    }
  
  fmpz_mod_poly_clear(red, fpctx);
  _fmpz_vec_clear(roots, nb);
  flint_free(mults);
  return res;
}
