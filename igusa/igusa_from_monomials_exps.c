
#include "igusa.h"

/* See igusa_base_exps */

/* Weight 20 */
static slong i_0020 = 0;
static slong i_1110 = 1;
static slong i_2001 = 2;
static slong i_5000 = 3;
static slong i_2200 = 4;

/* Weight 30 */ /*
static slong i_0500 = 0;
static slong i_3300 = 1;
static slong i_0030 = 2;
static slong i_5010 = 3;
static slong i_0301 = 4;
static slong i_1120 = 5;
static slong i_2011 = 6; */

/* Weight 60 */ /*
static slong i_0060 = 0;
static slong i_0005 = 1;
static slong i_0530 = 2;
static slong i_0204 = 3;
static slong i_3004 = 4;
static slong i_0X00 = 5;
static slong i_3800 = 6;
static slong i_5040 = 7;
static slong i_1150 = 8;
static slong i_3203 = 9; */

void igusa_from_monomials_exps(slong* e4, slong* e6, slong* e10, slong* e12,
			       int z4, int z6, int z10, int z12, slong wt)
{
  slong nb = igusa_nb_base_monomials(wt);
  slong k;

  for (k = 0; k < nb; k++)
    {
      e4[k] = 0;
      e6[k] = 0;
      e10[k] = 0;
      e12[k] = 0;
    }
  
  if (z10 && z12)
    {	      
      flint_printf("(igusa_rec_exps) Error: chi10, chi12 cannot be simultaneously zero\n");
      fflush(stdout);
      flint_abort();
    }
 
  if (wt == 20)
    {
      if (z4)
	{
	  flint_printf("(igusa_rec_exps) Error: psi4 cannot be zero in case of weight 20\n");
	  fflush(stdout);
	  flint_abort();	  
	}
      if (!z10)
	{
	  e4[i_0020] = 1;
	  e4[i_5000] = 1;

	  e6[i_0020] = 1;
	  e6[i_1110] = 1;
	  e6[i_5000] = 1;

	  e10[i_5000] = 2;
	  e10[i_0020] = 3;

	  e12[i_5000] = 2;
	  e12[i_2001] = 1;;
	  e12[i_0020] = 3;
	}
      else if (z10 && !z6)
	{
	  e4[i_5000] = 1;
	  e4[i_2200] = 1;
	  
	  e6[i_2200] = 2;
	  e6[i_5000] = 1;
	  
	  e12[i_2001] = 1;
	  e12[i_2200] = 3;
	  e12[i_5000] = 2;
	}
      else /* z10 and z6 */	    
	{
	  e4[i_5000] = 1;
	  
	  e12[i_5000] = 2;
	  e12[i_2001] = 1;
	}
    }

  else if (wt == 30)
    {

    }

  else if (wt == 60)
    {

    }
  
  else
    {
      flint_printf("(igusa_from_monomials) Error: weight %wd not implemented\n", wt);
      fflush(stdout);
      flint_abort();
    }
}
