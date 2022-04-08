
#include "modular.h"

/* field is assumed to be defined by a polynomial dividing num1, so that x is a root */
int modeq_isog_invariants_nf(fmpq_poly_struct* j, const fmpz_poly_struct* num_vec,
			     slong nb, const fmpz_poly_t field)
{
  fmpz_poly_t deriv;
  fmpq_poly_t mod_fmpq;
  fmpq_poly_t G;
  fmpq_poly_t S, T;
  fmpq_poly_t inverse;
  fmpq_t coeff;
  fmpq_poly_t v;
  slong k;
  int res;

  fmpz_poly_init(deriv);
  fmpq_poly_init(G);
  fmpq_poly_init(S);
  fmpq_poly_init(T);
  fmpq_poly_init(mod_fmpq);
  fmpq_poly_init(inverse);
  fmpq_init(coeff);
  fmpq_poly_init(v);
  
  fmpq_poly_zero(&j[0]);
  fmpq_poly_set_coeff_si(&j[0], 1, 1);
  fmpz_poly_derivative(deriv, &num_vec[0]);
  fmpq_poly_set_fmpz_poly(mod_fmpq, field);
  fmpq_poly_set_fmpz_poly(inverse, deriv);
  fmpq_poly_xgcd(G, S, T, inverse, mod_fmpq); /* Inverse of deriv should be S/G */
  
  if (fmpq_poly_degree(G) != 0) res = 0;
  else
    {
      res = 1;
      fmpq_poly_get_coeff_fmpq(coeff, G, 0);
      fmpq_poly_scalar_div_fmpq(inverse, S, coeff); 
      for (k = 1; k < nb; k++)
	{
	  /* Reduce num_vec[k] in number field */
	  fmpq_poly_set_fmpz_poly(v, &num_vec[k]);
	  fmpq_poly_rem(v, v, mod_fmpq);
	  fmpq_poly_mul(v, v, inverse);
	  fmpq_poly_rem(v, v, mod_fmpq);
	  fmpq_poly_set(&j[k], v);
	} 
    }

  fmpz_poly_clear(deriv);
  fmpq_poly_clear(G);
  fmpq_poly_clear(S);
  fmpq_poly_clear(T);
  fmpq_poly_clear(mod_fmpq);
  fmpq_poly_clear(inverse);
  fmpq_clear(coeff);
  fmpq_poly_clear(v);

  return res;
}
