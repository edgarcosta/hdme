
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_coset....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 3 * arb_test_multiplier(); iter++)
    {
      slong bits = 4;
      slong ell = n_randprime(state, bits, 1); /* Less than 2^bits */
      slong n = siegel_nb_cosets(ell);
      sp2gz_struct mats[4096]; /* Length is 2^(3*bits) */
      fmpz_mat_t temp1, temp2;
      fmpz_t det;
      slong k, i;

      fmpz_mat_init(temp1, 4, 4);
      fmpz_mat_init(temp2, 4, 4);
      fmpz_init(det);
      for (k = 0; k < n; k++) sp2gz_init(&mats[k], 2);
      
      for (k = 0; k < n; k++) siegel_coset(&mats[k], k, ell);
      /* Check:
	 - Determinant ell^2
	 - Lower left is zero
	 - All distinct */
      /* The matrices we have are of the form:
	 etaR * [1, 0; 0, l] * eta
	 We want to know there is no symplectic matrix such that
	 g * ... = ...
	 i.e., etaR * [100l] * eta * eta'^{-1} * [l001] * etaR^{-1} is NOT divisible by l. */
      for (k = 0; k < n; k++)
	{
	  if (!fmpz_mat_is_zero(&(&mats[k])->c))
	    {
	      flint_printf("FAIL (nonzero c)\n");
	      flint_printf("ell = %wd, k = %wd\n", ell, k);
	      sp2gz_print(&mats[k]);
	      fflush(stdout);
	      flint_abort();
	    }
	  sp2gz_get_mat(temp1, &mats[k]);
	  fmpz_mat_det(det, temp1);
	  if (!fmpz_equal_si(det, ell*ell))
	    {
	      flint_printf("FAIL (determinant)\n");
	      flint_printf("ell = %wd, k = %wd\n", ell, k);
	      sp2gz_print(&mats[k]);
	      fflush(stdout);
	      flint_abort();
	    }
	}
      for (k = 0; k < n; k++)
	{
	  for (i = k+1; i < n; i++)
	    {
	      sp2gz_get_mat(temp1, &mats[i]);
	      fmpz_mat_inv(temp1, det, temp1); /* det is the denominator */
	      if (!fmpz_equal_si(det, ell)) /* Then it's l^2: divide by l */
		{
		  fmpz_mat_scalar_divexact_si(temp1, temp1, ell);
		}
	      sp2gz_get_mat(temp2, &mats[k]);
	      fmpz_mat_mul(temp2, temp2, temp1);
	      fmpz_set_si(det, ell);
	      fmpz_mat_scalar_smod(temp2, temp2, det);
	      
	      if (fmpz_mat_is_zero(temp2))
		{
		  flint_printf("FAIL (not distinct)\n");
		  flint_printf("ell = %wd, k = %wd, i = %wd\n", ell, k, i);
		  fmpz_mat_print_pretty(temp1);
		  sp2gz_print(&mats[k]);
		  sp2gz_print(&mats[i]);
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}
	 
      fmpz_mat_clear(temp1);
      fmpz_mat_clear(temp2);
      fmpz_clear(det);
      for (k = 0; k < n; k++) sp2gz_clear(&mats[k]);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
