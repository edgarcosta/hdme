
#include "hecke.h"

void hecke_eigenvalue_eisenstein_p(fmpz_t eig, slong k, slong p)
{
  fmpz_t pp;
  fmpz_t temp;

  fmpz_init(pp);
  fmpz_init(temp);
  
  if (k < 4)
    {
      flint_printf("(hecke_eigenvalue_eisenstein_p) No Eisenstein series of weight %wd\n", k);
      fflush(stdout);
      flint_abort();
    }

  fmpz_set_si(eig, 1);
  fmpz_set_si(pp, p);
  
  fmpz_pow_ui(temp, pp, k-2);
  fmpz_add(eig, eig, temp);
  fmpz_pow_ui(temp, pp, k-2);
  fmpz_add(eig, eig, temp);
  fmpz_pow_ui(temp, pp, 2*k-3);
  fmpz_add(eig, eig, temp);

  fmpz_clear(pp);
  fmpz_clear(temp);
}
