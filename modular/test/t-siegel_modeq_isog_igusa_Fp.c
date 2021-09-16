
#include <string.h>
#include "modular.h"

int main()
{
  slong iter;
  
  flint_printf("siegel_modeq_isog_igusa_Fp....");
  fflush(stdout);

  for (iter = 0; iter < 4; iter++) /* 4 -> 6 to test with ell=7 */
    {
      
      fmpz_t p;
      slong ell;
      fmpz_poly_t crv1, crv2;
      fmpz_t coeff;
      slong crv1_coeffs[7], crv2_coeffs[7];
      fmpz* j1;
      fmpz* j2;
      fmpq* j1_Q;
      fmpq* j2_Q;
      slong nb_roots = 0;
      slong max_len;
      slong k;
      
      fmpz_mod_ctx_t ctx;
      fmpz_mod_poly_t pol1, pol2, pol3;
      fmpz* roots;
      slong* mults;
      int res = 1;

      /* Init things independent of iter */
      fmpz_init(p);
      fmpz_poly_init(crv1);
      fmpz_poly_init(crv2);
      fmpz_init(coeff);
      j1 = _fmpz_vec_init(3);
      j2 = _fmpz_vec_init(3);
      j1_Q = _fmpq_vec_init(3);
      j2_Q = _fmpq_vec_init(3);

      /* Set ell, p, curve coefficients */
      if (iter == 0)
	{
	  slong c1[7] = {2998, 1069, 4063, 4297, 4278, 272, 0};
	  slong c2[7] = {655, 1453, 4588, 3115, 2322, 695, 0};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 3;
	  fmpz_set_si(p, 5261);
	}	  
      if (iter == 1)
	{
	  slong c1[7] = {608155867,
			 2233670825,
			 671225738,
			 1092030439,
			 2028583210,
			 48872812,
			 1774507961};
	  slong c2[7] = {825611194,
			 183139891,
			 2422852850,
			 87910848,
			 1720123333,
			 2286039407,
			 1927466494};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 3;
	  fmpz_set_si(p, 2534267893);
	}
      if (iter == 2) 
	{
	  continue; /* Double root? */
	  /*
	  slong c1[7] = {14, 5, 59, 91, 71, 27, 0};
	  slong c2[7] = {51, 7, 20, 38, 26, 29, 0};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 5;
	  fmpz_set_si(p, 101); */
	}
      if (iter == 3)
	{
	  slong c1[7] = {2690010753,
			 159741347,
			 2805510580,
			 2536585478,
			 3653091983,
			 2420332800, 0};
	  slong c2[7] = {1266950302,
			 3172892980,
			 1209100179,
			 3748957676,
			 2616936853,
			 4076826784, 0};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 5;
	  fmpz_set_si(p, 4294967311);
	}
      if (iter == 4)
	{
	  slong c1[7] = {1308, 7948, 5411, 2876, 471, 4826, 0};
	  slong c2[7] = {4087, 1845, 7103, 7011, 7699, 7218, 0};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 7;
	  fmpz_set_si(p, 10009);
	}
      if (iter == 5)
	{
	  slong c1[7] = {3241061457,
			 1342046779,
			 448279016,
			 471351782,
			 1698662093,
			 393356368, 0};
	  slong c2[7] = {474327543,
			 306153493,
			 580698082,
			 2005208933,
			 2231412358,
			 2171506943, 0};
	  memcpy(crv1_coeffs, c1, 7 * sizeof(slong));
	  memcpy(crv2_coeffs, c2, 7 * sizeof(slong));
	  ell = 7;
	  fmpz_set_si(p, 3452678353);
	}
      /* Init things depending on iter */
      fmpz_mod_ctx_init(ctx, p);
      fmpz_mod_poly_init(pol1, ctx);
      fmpz_mod_poly_init(pol2, ctx);
      fmpz_mod_poly_init(pol3, ctx);
      
      max_len = siegel_nb_cosets(ell);
      roots = _fmpz_vec_init(max_len);
      mults = flint_malloc(max_len * sizeof(slong));

      /* Set curves */
      fmpz_poly_zero(crv1);
      fmpz_poly_zero(crv2);
      for (k = 0; k <= 6; k++)
	{
	  fmpz_set_si(coeff, crv1_coeffs[k]);
	  fmpz_poly_set_coeff_fmpz(crv1, k, coeff);
	  fmpz_set_si(coeff, crv2_coeffs[k]);
	  fmpz_poly_set_coeff_fmpz(crv2, k, coeff);	  
	}

      igusa_from_curve_fmpz(j1_Q, crv1);
      igusa_from_curve_fmpz(j2_Q, crv2);
      for (k = 0; k < 3; k++)
	{
	  fmpq_mod_fmpz(&j1[k], &j1_Q[k], p);
	  fmpq_mod_fmpz(&j2[k], &j2_Q[k], p);
	}

      /* Check roots of modular equations */
      siegel_modeq_eval_Fp(pol1, pol2, pol3, j1, ell, ctx);
      siegel_modeq_roots_Fp(&nb_roots, roots, mults, pol1, ctx);
      res = 0;
      for (k = 0; k < nb_roots; k++)
	{
	  if (fmpz_equal(&roots[k], &j2[0])) res = 1;
	}
      if (!res)
	{
	  flint_printf("FAIL (root not found)\n");
	  fmpz_print(&j2[0]); flint_printf("\n");
	  flint_printf("iter = %wd; nb_roots = %wd\n", iter, nb_roots);
	  for (k = 0; k < nb_roots; k++)
	    {
	      fmpz_print(&roots[k]); flint_printf("\n");
	    }
	  fmpz_poly_print_pretty(crv1, "x"); flint_printf("\n");
	  fmpz_poly_print_pretty(crv2, "x"); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      /* Check other Igusa invariants, using j1 as temp */
      res = siegel_modeq_isog_igusa_Fp(j1, pol1, pol2, pol3, &j2[0], ctx);
      if (!res)
	{
	  flint_printf("FAIL (isog_igusa)\n");
	  fmpz_print(&j2[0]);
	  flint_printf("iter = %wd\n", iter);
	  fflush(stdout);
	  flint_abort();
	}
      for (k = 0; k < 3; k++)
	{
	  if (!fmpz_equal(&j1[k], &j2[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (igusa)\n");
	  flint_printf("iter = %wd\n", iter);
	  for (k = 0; k < 3; k++)
	    {
	      fmpz_print(&j1[k]); flint_printf("\n");
	      fmpz_print(&j2[k]); flint_printf("\n");
	    }
	  fmpz_poly_print_pretty(crv1, "x"); flint_printf("\n");
	  fmpz_poly_print_pretty(crv2, "x"); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      /* Clear all */
      fmpz_clear(p);
      fmpz_poly_clear(crv1);
      fmpz_poly_clear(crv2);
      fmpz_clear(coeff);
      _fmpz_vec_clear(j1, 3);
      _fmpz_vec_clear(j2, 3);
      _fmpq_vec_clear(j1_Q, 3);
      _fmpq_vec_clear(j2_Q, 3);
      
      fmpz_mod_poly_clear(pol1, ctx);
      fmpz_mod_poly_clear(pol2, ctx);
      fmpz_mod_poly_clear(pol3, ctx);
      fmpz_mod_ctx_clear(ctx);
      _fmpz_vec_clear(roots, max_len);
      flint_free(mults);
    }
  
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

