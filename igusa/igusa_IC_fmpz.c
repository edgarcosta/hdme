
#include "igusa.h"

void igusa_IC_fmpz(fmpz* IC, fmpz* I)
{
  fmpz_t g, I2, I6;
  fmpz* resc;
  int r;

  fmpz_init(g);
  fmpz_init(I2);
  fmpz_init(I6);
  resc = _fmpz_vec_init(4);

  fmpz_gcd(g, cov_I12(I), cov_I10(I));
  fmpz_divexact(g, cov_I10(I), g);
  cov_rescale_fmpz(resc, I, g); /* Ensure I2 is an integer */
  
  r = igusa_I2_fmpz(I2, resc);
  if (!r)
    {
      flint_printf("(igusa_IC_fmpz) Error: I2 is not integral\n");
      fflush(stdout);
      flint_abort();
    }

  r = igusa_I6_fmpz(I6, resc);
  if (!r)
    {
      cov_rescale_fmpz_si(resc, resc, 3);
      igusa_I2_fmpz(I2, resc);
      igusa_I6_fmpz(I6, resc);
    }

  fmpz_set(&IC[0], I2);
  fmpz_set(&IC[1], cov_I4(resc));
  fmpz_set(&IC[2], I6);
  fmpz_set(&IC[3], cov_I10(resc));
  /* No need to normalize, if I was previously normalized itself. */

  fmpz_clear(g);
  fmpz_clear(I2);
  fmpz_clear(I6);
  _fmpz_vec_clear(resc, 4);    
}
