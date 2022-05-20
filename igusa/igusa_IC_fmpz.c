
#include "igusa.h"

void igusa_IC_fmpz(fmpz* IC, fmpz* I)
{
  fmpz_t I2, I6;
  fmpz* resc;
  slong weights[4] = IGUSA_HALFWEIGHTS;
  int r;

  fmpz_init(I2);
  fmpz_init(I6);
  resc = _fmpz_vec_init(4);

  cov_rescale_fmpz(resc, I, igusa_I10(I), 4, weights); /* Ensure I2 is an integer */
  cov_rescale_fmpz_si(resc, resc, 3, 4, weights); /* Ensure I6 is an integer */
  
  r = igusa_I2_fmpz(I2, resc) && igusa_I6_fmpz(I6, resc);
  if (!r)
    {
      flint_printf("(igusa_IC_fmpz) Error: I2 or I6 is not integral\n");
      fflush(stdout);
      flint_abort();
    }

  fmpz_set(&IC[0], I2);
  fmpz_set(&IC[1], igusa_I4(resc));
  fmpz_set(&IC[2], I6);
  fmpz_set(&IC[3], igusa_I10(resc));
  cov_normalize_fmpz(IC, IC, 4, weights);

  fmpz_clear(I2);
  fmpz_clear(I6);
  _fmpz_vec_clear(resc, 4);    
}
