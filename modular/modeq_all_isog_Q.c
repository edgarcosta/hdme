
#include "modular.h"

void modeq_all_isog_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
		      slong* exp_array, const modeq_t E, const modeq_ctx_t ctx)
{
  slong nb = modeq_degree(E);
  slong* mults;
  fmpq* roots;
  slong k;

  mults = flint_malloc(nb * sizeof(slong));
  roots = _fmpq_vec_init(nb);
  
  pol_roots_Q(nb_roots, roots, mults, modeq_equation(E));
  *nb_M = modeq_ctx_nb(ctx);
  for (k = 0; k < *nb_M; k++)
    {
      cov_monomial_degrees(&exp_array[4*k], modeq_ctx_monomial(ctx, k),
			   modeq_ctx_ctx(ctx));
    }
  for (k = 0; k < *nb_roots; k++)
    {
      modeq_isog_monomials_Q(&all_M[(*nb_M) * k], E, &roots[k], mults[k]);
    }

  flint_free(mults);
  _fmpq_vec_clear(roots, nb);
}
  
