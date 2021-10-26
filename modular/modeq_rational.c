
#include "modular.h"

int modeq_rational(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* pol_vec_acb,
		   slong degree, slong nb, slong prec)
{
  fmpz* den_vec;
  fmpq_poly_struct* pol_vec;
  slong k;
  int success = 1;

  den_vec = flint_malloc(nb * sizeof(fmpz));
  pol_vec = flint_malloc(nb * sizeof(fmpq_poly_struct));
  for (k = 0; k < nb; k++)
    {
      fmpz_init(&den_vec[k]);
      fmpq_poly_init(&pol_vec[k]);
    }

  for (k = 0; k < nb; k++)
    {
      if (success)
	{
	  success = modeq_rational_poly(&pol_vec[k], &pol_vec_acb[k],
					degree, prec);
	}
    }

  if (success)
    { 
      for (k = 0; k < nb; k++)
	{
	  fmpq_poly_get_denominator(&den_vec[k], &pol_vec[k]);
	}
      _fmpz_vec_lcm(den, den_vec, nb);
      for (k = 0; k < nb; k++)
	{
	  fmpq_poly_scalar_mul_fmpz(&pol_vec[k], &pol_vec[k], den);
	  fmpq_poly_get_numerator(&num_vec[k], &pol_vec[k]);
	}
    }

  for (k = 0; k < nb; k++)
    {
      fmpz_clear(&den_vec[k]);
      fmpq_poly_clear(&pol_vec[k]);
    }
  flint_free(den_vec);
  flint_free(pol_vec);
  
  return success;
}
