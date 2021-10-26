
#include "modular.h"

int modeq_reduce(fmpz_mod_poly_struct* red_vec, const fmpz_poly_struct* num_vec,
		 const fmpz_t den, slong nb, const fmpz_mod_ctx_t ctx)
{
  int res;
  fmpz_t one;
  fmpz_t invden;
  int v = MODEQ_VERBOSE;
  slong k;
  
  fmpz_init(one);
  fmpz_init(invden);
  
  fmpz_one(one);
  res = fmpz_mod_divides(invden, one, den, ctx);
  
  if (v && !res) flint_printf("(modeq_reduce) Denominator vanishes\n");

  if (res)
    {
      for (k = 0; k < nb; k++)
	{
	  fmpz_mod_poly_set_fmpz_poly(&red_vec[k], &num_vec[k], ctx);
	  fmpz_mod_poly_scalar_mul_fmpz(&red_vec[k], &red_vec[k], invden, ctx);
	}
    }

  fmpz_clear(one);
  fmpz_clear(invden);
  return res;  
}
