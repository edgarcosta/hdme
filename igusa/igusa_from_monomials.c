
#include "igusa.h"

/* Semantic access: see igusa_base_exps */

/* Weight 20 */
static fmpz* mon_0020(fmpz* mon) {return &mon[0]};
static fmpz* mon_1110(fmpz* mon) {return &mon[1]};
static fmpz* mon_2001(fmpz* mon) {return &mon[2]};
static fmpz* mon_5000(fmpz* mon) {return &mon[3]};
static fmpz* mon_2200(fmpz* mon) {return &mon[4]};

/* Weight 30 */
static fmpz* mon_0500(fmpz* mon) {return &mon[0]};
static fmpz* mon_3300(fmpz* mon) {return &mon[1]};
static fmpz* mon_0030(fmpz* mon) {return &mon[2]};
static fmpz* mon_5010(fmpz* mon) {return &mon[3]};
static fmpz* mon_0301(fmpz* mon) {return &mon[4]};
static fmpz* mon_1120(fmpz* mon) {return &mon[5]};
static fmpz* mon_2011(fmpz* mon) {return &mon[6]};

/* Weight 60 */
static fmpz* mon_0060(fmpz* mon) {return &mon[0]};
static fmpz* mon_0005(fmpz* mon) {return &mon[1]};
static fmpz* mon_0530(fmpz* mon) {return &mon[2]};
static fmpz* mon_0204(fmpz* mon) {return &mon[3]};
static fmpz* mon_3004(fmpz* mon) {return &mon[4]};
static fmpz* mon_0X00(fmpz* mon) {return &mon[5]};
static fmpz* mon_3800(fmpz* mon) {return &mon[6]};
static fmpz* mon_5040(fmpz* mon) {return &mon[7]};
static fmpz* mon_1150(fmpz* mon) {return &mon[8]};
static fmpz* mon_3203(fmpz* mon) {return &mon[9]};

static void igusa_20(fmpz* I, fmpz* mon)
{
  fmpz_t m1, m2, m3, m4;
  fmpz_t temp;
  fmpz* res;
  slong wts[4] = IGUSA_HALFWEIGHTS;
  
  fmpz_init(m1);
  fmpz_init(m2);
  fmpz_init(m3);
  fmpz_init(m4);
  fmpz_init(temp);
  res = _fmpz_vec_init(4);

  _fmpz_vec_zero(res, 4);
  
  if (fmpz_is_zero(mon_5000(mon)))
    {
      flint_printf("(igusa_from_monomials) Error: psi4 cannot be zero in case of weight 20\n");
      fflush(stdout);
      flint_abort();
    }

  if (!fmpz_is_zero(mon_0020(mon)))
    {
      /* chi10 is nonzero; set chi10 = chi4^2 */
      fmpz_set(m1, mon_2001(mon));
      fmpz_set(m2, mon_0020(mon));
      fmpz_set(m3, mon_1110(mon));
      fmpz_set(m4, mon_5000(mon));

      fmpz_mul(igusa_psi4(res), m2, m4);
      fmpz_mul(igusa_psi6(res), m2, m3);
      fmpz_mul(igusa_psi6(res), igusa_psi6(res), m4);
      fmpz_mul(igusa_chi10(res), m4, m4);
      fmpz_pow_ui(temp, m2, 3);
      fmpz_mul(igusa_chi10(res), igusa_chi10(res), temp);
      fmpz_mul(igusa_chi12(res), m4, m4);
      fmpz_mul(igusa_chi12(res), igusa_chi12(res), m1);
      fmpz_pow_ui(temp, m2, 3);
      fmpz_mul(igusa_chi12(res), igusa_chi12(res), temp);
      
      cov_normalize_fmpz(res, res, 4, weights);
    }
  else
    {
      /* chi10 is zero, so chi12 is nonzero */
      if (fmpz_is_zero(mon_2001(mon)))
	{
	  flint_printf("(igusa_from_monomials) Error: chi10, chi12 are simultaneously zero\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      if (fmpz_is_zero(mon_2200(mon)))
	{
	  /* psi6, chi10 are zero: set chi12 = psi4^2, then rescale */
	  fmpz_set(m1, mon_5000(mon));
	  fmpz_set(m2, mon_2001(mon));
	  
	  fmpz_set(igusa_psi4(res), m1);
	  fmpz_mul(igusa_chi12(res), m1, m1);
	  fmpz_mul(igusa_chi12(res), igusa_chi12(res), m2);

	  wts[0] = 1;
	  wts[1] = 0;
	  wts[2] = 0;
	  wts[3] = 3;
	  cov_normalize_fmpz(res, res, 4, wts);
	}
      else
	{
	  /* Set psi6 = psi4, then rescale */
	  fmpz_set(m1, mon_5000(mon));
	  fmpz_set(m2, mon_2001(mon));
	  fmpz_set(m3, mon_2200(mon));
	  
	  fmpz_mul(igusa_psi4(res), m1, m3);
	  fmpz_mul(igusa_psi6(res), m3, m3);
	  fmpz_mul(igusa_psi6(res), igusa_psi6(res), m1);
	  fmpz_pow_ui(igusa_chi12(res), m3, 3);
	  fmpz_mul(temp, m1, m1);
	  fmpz_mul(temp, temp, m2);
	  fmpz_mul(igusa_chi12(res), igusa_chi12(res), temp);
	  
	  cov_normalize_fmpz(res, res, 4, wts);
	}
    }

  _fmpz_vec_set(I, res, 4);

  fmpz_clear(m1);
  fmpz_clear(m2);
  fmpz_clear(m3);
  fmpz_clear(m4);
  fmpz_clear(temp);
  _fmpz_vec_clear(res, 4);
}

static void igusa_30(fmpz* I, fmpz* mon)
{

}

static void igusa_60(fmpz* I, fmpz* mon)
{

}


void igusa_from_monomials(fmpz* I, fmpz* mon, slong wt)
{
  if (wt == 20) igusa_20(I, mon);
  else if (wt == 30) igusa_30(I, mon);
  else if (wt == 60) igusa_60(I, mon);
  else
    {
      flint_printf("(igusa_from_monomials) Error: weight %wd not implemented\n", wt);
      fflush(stdout);
      flint_abort();
    }
}

