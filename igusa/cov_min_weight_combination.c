
#include "igusa.h"

void cov_min_weight_combination(slong* wt, slong* i1, slong* i2,
				slong* e1, slong* e2, fmpz* I)
{
  *wt = 2; /* Default */
  if (!fmpz_is_zero(cov_I4(I)))
    {      
      *i1 = 0;
      if (!fmpz_is_zero(cov_I6prime(I)))
	{
	  *i2 = 1;
	  *e1 = -1;
	  *e2 = 1;
	}
      else if (!fmpz_is_zero(cov_I10(I)))
	{
	  *i2 = 2;
	  *e1 = -2;
	  *e2 = 1;
	}
      else
	{
	  *wt = 4;
	  *i2 = *i1;
	  *e2 = 0;
	}
    }
  else if (!fmpz_is_zero(cov_I6prime(I)))
    {
      *i1 = 1;
      if (!fmpz_is_zero(cov_I10(I)))
	{
	  *i2 = 2;
	  *e1 = 2;
	  *e2 = -1;
	}
      else
	{
	  *wt = 6;
	  *i2 = *i1;
	  *e1 = 1;
	  *e2 = 0;
	}
    }
  else if (!fmpz_is_zero(cov_I10(I)))
    {
      *i1 = 2;
      if (!fmpz_is_zero(cov_I12(I)))
	{
	  *i2 = 3;
	  *e1 = -1;
	  *e2 = 1;
	}
      else
	{
	  *wt = 10;
	  *i2 = *i1;
	  *e1 = 1;
	  *e2 = 0;
	}
    }
  else if (!fmpz_is_zero(cov_I12(I)))
    {
      *wt = 12;
      *i1 = 3;
      *i2 = *i1;
      *e1 = 1;
      *e2 = 0;
    }
  else
    {
      flint_printf("(cov_min_weight_combination) Error: all entries are zero\n");
      fflush(stdout);
      flint_abort();
    }
}
