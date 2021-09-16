
#include "modular.h"

void siegel_coset(sp2gz_t m, slong k, slong ell)
{
  slong a, b, c;
  sp2gz_t eta, etaR;
  fmpz_mat_t temp;
  slong u, v, w, x, y, z;
  
  sp2gz_init(eta, 2);
  sp2gz_init(etaR, 2);
  fmpz_mat_init(temp, 4, 4);
  
  if ((k < 0) || (k >= siegel_nb_cosets(ell)))
    {
      flint_printf("(siegel_coset) No matrix numbered %wd\n", k);
      fflush(stdout);
      flint_abort();
    }
  
  if (k < n_pow(ell, 3))
    {
      /* Case T_1(a, b, c) */
      a = k % ell;
      k = k / ell;
      b = k % ell;
      k = k / ell;
      c = k % ell;
      
      fmpz_mat_one(temp);
      fmpz_mat_neg(temp, temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 2), a);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 3), c);
      sp2gz_set_mat(eta, temp);
      
      sp2gz_one(etaR);
    }
  else if (k < n_pow(ell, 3) + ell + (ell - 1) * ell)
    {
      /* Case T_2(a, b, c) */
      if (k < n_pow(ell, 3) + ell)
	{
	  a = 0;
	  b = 0;
	  c = k - n_pow(ell, 3);

	  /* Special etaR in this case */
	  if (c == 0)
	    {
	      sp2gz_J(etaR);
	    }
	  else
	    {
	      n_xgcd((ulong*) &v, (ulong*) &u, ell, c); 
	      u = -u; /* uc + lv = 1 */
	      fmpz_mat_zero(temp);
	      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), u);
	      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), v);
	      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), 1);
	      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), -ell);
	      fmpz_set_si(fmpz_mat_entry(temp, 2, 3), c);
	      fmpz_set_si(fmpz_mat_entry(temp, 3, 0), -1);
	      sp2gz_set_mat(etaR, temp);
	    }
	}
      else
	{
	  a = (k - n_pow(ell, 3) - ell) % (ell-1);
	  a += 1;
	  b = (k - n_pow(ell, 3) - ell) % ell;
	  c = (b*b) % ell;
	  c *= n_invmod(a, ell);
	  c = c % ell;

	  n_xgcd((ulong*) &v, (ulong*) &u, ell, a);
	  u = -u; /* au + lv = 1 */
	  w = - b*u;
	  w = w % ell; /* w = -b/a mod ell, u = 1/a mod ell */
	  x = - (b + w*a) / ell;
	  y = - (c + w*b) / ell;
	  z = u*x - w*v;

	  fmpz_mat_zero(temp);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 0), u);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 2), v);
	  fmpz_set_si(fmpz_mat_entry(temp, 0, 3), z);
	  fmpz_set_si(fmpz_mat_entry(temp, 1, 3), -1);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -ell);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 2), a);
	  fmpz_set_si(fmpz_mat_entry(temp, 2, 3), b);
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 0), w);
	  fmpz_one(fmpz_mat_entry(temp, 3, 1));
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 2), x);
	  fmpz_set_si(fmpz_mat_entry(temp, 3, 3), y);
	  sp2gz_set_mat(etaR, temp);
	}
      /* eta is always given by the same formula */
      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 0), -a);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), -b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 0), -b);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 1), -c);
      fmpz_one(fmpz_mat_entry(temp, 0, 2));
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      sp2gz_set_mat(eta, temp);
    }
  else if (k < n_pow(ell, 3) + ell*ell + ell)
    {
      /* Case T_3(a) */
      a = k - n_pow(ell, 3) - ell*ell;

      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 0), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 1), -a);
      fmpz_set_si(fmpz_mat_entry(temp, 1, 2), -a);
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 2), -1);
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      sp2gz_set_mat(eta, temp);

      fmpz_mat_zero(temp);
      fmpz_set_si(fmpz_mat_entry(temp, 0, 3), -1);
      fmpz_one(fmpz_mat_entry(temp, 1, 0));
      fmpz_one(fmpz_mat_entry(temp, 2, 1));
      fmpz_one(fmpz_mat_entry(temp, 3, 2));
      sp2gz_set_mat(etaR, temp);
    }
  else
    {
      /* Case T_4 */
      fmpz_mat_zero(temp);
      fmpz_one(fmpz_mat_entry(temp, 0, 1));
      fmpz_one(fmpz_mat_entry(temp, 1, 1));
      fmpz_one(fmpz_mat_entry(temp, 1, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), -1);
      fmpz_one(fmpz_mat_entry(temp, 2, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 3, 0), -1);
      fmpz_one(fmpz_mat_entry(temp, 3, 1));
      sp2gz_set_mat(eta, temp);

      fmpz_mat_zero(temp);
      fmpz_one(fmpz_mat_entry(temp, 0, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 1, 1), ell);
      fmpz_one(fmpz_mat_entry(temp, 1, 3));
      fmpz_set_si(fmpz_mat_entry(temp, 2, 0), -ell);
      fmpz_set_si(fmpz_mat_entry(temp, 2, 1), ell);
      fmpz_one(fmpz_mat_entry(temp, 2, 2));
      fmpz_one(fmpz_mat_entry(temp, 2, 3));
      fmpz_one(fmpz_mat_entry(temp, 3, 0));
      fmpz_set_si(fmpz_mat_entry(temp, 3, 1), -1);
      sp2gz_set_mat(etaR, temp);
    }
  
  if (!sp2gz_is_correct(eta) || !sp2gz_is_correct(etaR))
    {
      sp2gz_print(eta);
      sp2gz_print(etaR);
      flint_printf("k = %wd, a = %wd, b = %wd, c = %wd, ell = %wd\n", k, a, b, c, ell);
      flint_printf("eta is correct? %wd; etaR is correct? %wd\n",
		   sp2gz_is_correct(eta), sp2gz_is_correct(etaR));
      fflush(stdout);
      flint_abort();
    }
  /* Compute transformation directly: warning, they are only general
     symplectic */
  fmpz_mat_scalar_mul_si(&eta->c, &eta->c, ell);
  fmpz_mat_scalar_mul_si(&eta->d, &eta->d, ell);
  sp2gz_mul(m, etaR, eta);
  
  sp2gz_clear(eta);
  sp2gz_clear(etaR);
  fmpz_mat_clear(temp);
}
