
#include "hilbert.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_parametrize....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      slong delta;
      fmpq* rs;
      acb_t r_acb, s_acb;
      acb_ptr I;
      
      slong rs_bits = 5 + n_randint(state, 10);
      slong prec = 100 + n_randint(state, 1000);
      slong delta_max = 18;
      slong k;

      rs = _fmpq_vec_init(2);
      acb_init(r_acb);
      acb_init(s_acb);
      I = _acb_vec_init(4);
	
      for (delta = 5; delta < delta_max; delta++)
	{
	  if (hilbert_is_fundamental(delta))
	    {
	      for (k = 0; k < 2; k++)
		{
		  fmpq_randbits(&rs[k], state, rs_bits);
		}
	      acb_set_fmpq(r_acb, &rs[0], prec);
	      acb_set_fmpq(s_acb, &rs[1], prec);
	      hilbert_parametrize(I, r_acb, s_acb, delta, prec);
	    }
	}

      _fmpq_vec_clear(rs, 2);
      acb_clear(r_acb);
      acb_clear(s_acb);
      _acb_vec_clear(I, 4);      
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
