
#include "igusa.h"

/* Cf. Igusa, "On Siegel modular forms of genus two", p. 181 
   In the static functions, I is the Streng invariants. */

static void igusa_y(acb_t y1, acb_t y2, acb_srcptr I, slong prec)
{
  acb_pow_si(y1, &I[0], 3, prec);
  acb_div(y1, y1, &I[3], prec);
  acb_mul_si(y1, y1, n_pow(2,11) * 3, prec);

  acb_pow_si(y2, &I[1], 2, prec);
  acb_div(y2, y2, &I[3], prec);
  acb_mul_si(y2, y2, n_pow(2,13) * 3, prec);
}

static void igusa_ec_j1j2_gen(acb_ptr j, const acb_t y1, const acb_t y2, slong prec)
{
  acb_poly_t pol;
  acb_t c;

  acb_poly_init(pol);
  acb_init(c);
  
  acb_poly_set_coeff_si(pol, 2, 1);
  acb_poly_set_coeff_acb(pol, 0, y1);

  acb_sub(c, y2, y1, prec);
  acb_sub_si(c, c, n_pow(2,12) * n_pow(3,6), prec);
  acb_div_si(c, c, n_pow(2,6) * n_pow(3,3), prec);
  acb_poly_set_coeff_acb(pol, 1, c);

  /* Assume that thomae_roots will not throw for deg 2 polynomials? */
  thomae_roots(j, pol, prec);
  
  acb_poly_clear(pol);
  acb_clear(c);
}

static void igusa_ec_j1j2_zeroI4(acb_ptr j, const acb_t y2, slong prec)
{
  acb_zero(&j[0]);
  acb_div_si(&j[1], y2, -n_pow(2,6)*n_pow(3,3), prec); /* aka 1728 */
  acb_add_si(&j[1], &j[1], n_pow(2,6)*n_pow(3,3), prec);
}

static void igusa_ec_j1j2_zeroI6prime(acb_ptr j, const acb_t y1, slong prec)
{
  acb_set_si(&j[0], 1728);
  acb_div_si(&j[1], y1, 1728, prec);
}

static void igusa_ec_j1j2_zeroI4I6prime(acb_ptr j)
{
  acb_set_si(&j[0], 0);
  acb_set_si(&j[1], 1728);
}


void igusa_ec_j1j2(acb_ptr j, fmpz* I, slong prec)
{
  acb_ptr S_acb;
  slong k;
  acb_t y1, y2;
  fmpz* S;

  S_acb = _acb_vec_init(4);
  acb_init(y1);
  acb_init(y2);
  S = _fmpz_vec_init(4);

  igusa_streng_fmpz(S, I);
  
  for (k = 0; k < 4; k++) acb_set_fmpz(&S_acb[k], &S[k]);
  igusa_y(y1, y2, S_acb, prec);

  if (fmpz_is_zero(&S[0]) && fmpz_is_zero(&S[1]))
    {
      igusa_ec_j1j2_zeroI4I6prime(j);
    }
  else if (fmpz_is_zero(&S[0]))
    {
      igusa_ec_j1j2_zeroI4(j, y2, prec);
    }
  else if (fmpz_is_zero(&S[1]))
    {
      igusa_ec_j1j2_zeroI6prime(j, y1, prec);
    }
  else
    {
      igusa_ec_j1j2_gen(j, y1, y2, prec);
    }

  _acb_vec_clear(S_acb, 4);
  acb_clear(y1);
  acb_clear(y2);
  _fmpz_vec_clear(S, 4);
}

