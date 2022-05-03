
#include "modular.h"

int siegel_modeq_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j,
				   fmpq* j, slong ell)
{
  
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  slong* mults;
  fmpq* roots;
  slong nmax = siegel_nb_cosets(ell);
  slong k;
  int res;
  
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  roots = _fmpq_vec_init(nmax);
  mults = flint_malloc(nmax * sizeof(slong));
  
  res = siegel_modeq_eval_Q(num_vec, den, j, ell);
  
  if (res) modeq_roots_Q(nb_roots, roots, mults, &num_vec[0]);
      
  for (k = 0; (k < *nb_roots) && res; k++)
    {
      res = modeq_isog_invariants_Q(&all_isog_j[3*k], num_vec, &roots[k], 3);
    }
  
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  _fmpq_vec_clear(roots, nmax);
  flint_free(mults);
  return res;
}
