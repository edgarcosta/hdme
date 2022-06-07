
#include "hecke.h"

int main()
{
  slong iter;
  
  flint_printf("hecke_make_integral....");
  fflush(stdout);

  for (iter = 0; iter < 1; iter++)
    {
      fmpz* I;
      slong ell = 2;
      hecke_t H;
      slong nb = siegel_nb_T1_cosets(ell);
      acb_t psi4;
      slong w = 4;
      acb_t scal;
      acb_t c;
      slong k;
      int print = 1;
      slong prec = 1000;
      slong weights[4] = IGUSA_WEIGHTS;

      I = _fmpz_vec_init(4);
      ell = 2;
      hecke_init(H, nb);
      acb_init(psi4);
      acb_init(scal);
      acb_init(c);
            
      fmpz_set_si(&I[0], 108);
      fmpz_set_si(&I[1], 57);
      fmpz_set_si(&I[2], 2259);
      fmpz_set_si(&I[3], -31872);
      igusa_from_IC_fmpz(I, I);
      hecke_set_I_fmpz(H, I, prec);
      hecke_collect_T1(H, ell, prec);

      cov_find_rescaling(scal, hecke_I_tau(H), I, 4, weights, prec);

      for (k = 0; k < nb; k++)
	{
	  acb_set(psi4, igusa_psi4(hecke_I(H, k)));
	  acb_pow_si(c, scal, w, prec);
	  acb_div(psi4, psi4, c, prec);
	  acb_inv(c, hecke_stardet(H, k), prec);
	  acb_mul_si(c, c, n_pow(ell, 3), prec);
	  acb_pow_si(c, c, w, prec);
	  acb_mul(psi4, psi4, c, prec);
	  acb_mul_si(psi4, psi4, 6, prec);

	  if (print)
	    {
	      flint_printf("Is this integral?\n");
	      acb_printd(psi4, 20); flint_printf("\n");
	    }
	}

      _fmpz_vec_clear(I, 4);
      hecke_clear(H);
      acb_clear(psi4);
      acb_clear(scal);
      acb_clear(c);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
