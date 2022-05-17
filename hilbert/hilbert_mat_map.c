
#include "hilbert.h"

void hilbert_mat_map(fmpz_mat_t eta, const fmpz_poly_mat_t m, slong delta)
{
  slong prec = 100;
  int success = 0;

  acb_mat_t Ri, Rt;
  acb_mat_t M, N, P;
  slong j, k;

  acb_mat_init(Ri, 2, 2);
  acb_mat_init(Rt, 2, 2);
  acb_mat_init(M, 4, 4);
  acb_mat_init(N, 4, 4);
  acb_mat_init(P, 4, 4);

  while (!success)
    {
      success = 1;
      hilbert_R(Ri, delta, prec);
      acb_mat_transpose(Rt, Ri);
      acb_mat_inv(Ri, Ri, prec);
      
      acb_mat_zero(M);
      acb_mat_set_window(M, 0, 0, Rt);
      acb_mat_set_window(M, 2, 2, Ri);
      acb_mat_inv(P, M, prec);

      hilbert_sigma1(acb_mat_entry(N, 0, 0),
		     fmpz_poly_mat_entry(m, 0, 0), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 1, 1),
		     fmpz_poly_mat_entry(m, 0, 0), delta, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 0, 2),
		     fmpz_poly_mat_entry(m, 0, 1), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 1, 3),
		     fmpz_poly_mat_entry(m, 0, 1), delta, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 2, 0),
		     fmpz_poly_mat_entry(m, 1, 0), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 2, 1),
		     fmpz_poly_mat_entry(m, 1, 0), delta, prec);
      
      hilbert_sigma1(acb_mat_entry(N, 2, 2),
		     fmpz_poly_mat_entry(m, 1, 1), delta, prec);
      hilbert_sigma2(acb_mat_entry(N, 3, 3),
		     fmpz_poly_mat_entry(m, 1, 1), delta, prec);

      acb_mat_mul(M, M, N, prec);
      acb_mat_mul(M, M, P, prec);
      for (j = 0; j < 4; j++)
	{
	  for (k = 0; k < 4; k++)
	    {
	      success =
		success && acb_round(fmpz_mat_entry(eta, j, k),
				     acb_mat_entry(M, j, k));
	    }
	}
      prec *= 2;
      if (!success)
	{
	  flint_printf("(hilbert_mat_map) Info: increasing precision\n");
	  acb_mat_printd(M, 10); flint_printf("\n");
	}
    }
    
  acb_mat_clear(Ri);
  acb_mat_clear(Rt);
  acb_mat_clear(M);
  acb_mat_clear(N);
  acb_mat_clear(P);  
}
