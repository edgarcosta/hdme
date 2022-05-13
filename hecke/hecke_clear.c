
#include "hecke.h"

void hecke_clear(hecke_t H)
{
  slong k;
  
  acb_mat_clear(hecke_tau(H));
  _acb_vec_clear(hecke_t1t2(H), 2);
  fmpz_poly_clear(hecke_beta(H));
  
  for (k = 0; k < nb; k++) fmpz_mat_clear(hecke_coset(H, k));
  for (k = 0; k < nb; k++) acb_mat_clear(hecke_isog(H, k));
  for (k = 0; k < nb; k++) acb_mat_clear(hecke_isog(H, k));  

  flint_free(H->cosets);
  flint_free(H->isog);
  flint_free(H->stars);
  _acb_vec_clear(H->stardets, nb);
  _acb_vec_clear(H->theta2, 16*nb);
  _acb_vec_clear(H->I, 4*nb);  
}
