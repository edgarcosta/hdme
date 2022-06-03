
#include "hecke.h"

int main()
{
  slong wt;
  
  flint_printf("hecke_charpoly....");
  fflush(stdout);

  for (wt = 4; wt < 100; wt += 2)
    {      
      slong ell = 2;
      fmpz_poly_t pol;

      fmpz_poly_init(pol);
      hecke_charpoly(pol, ell, wt);

      flint_printf("Charpoly for weight %wd:\n", wt);
      fmpz_poly_print(pol); flint_printf("\n");

      fmpz_poly_clear(pol);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


