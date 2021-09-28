
#include "hilbert.h"

int main()
{

  slong delta;
  char** vars;
  int print = 1;

  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_t test;
  
  flint_printf("humbert_get_mpoly....");
  fflush(stdout);

  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(test, ctx);
  vars = humbert_vars_init();

  for (delta = 0; delta < 100; delta++)
    {
      if (hilbert_is_fundamental(delta) && (delta != 33))
	{
	  humbert_vars_set(vars, delta);
	  
	  humbert_get_mpoly(test, (const char**) vars, "A1", delta, ctx);
	  if (print)
	    {
	      flint_printf("A1(%wd) = ", delta);
	      fmpq_mpoly_print_pretty(test, (const char**) vars, ctx);
	      flint_printf("\n");
	    }
	    
	  humbert_get_mpoly(test, (const char**) vars, "A", delta, ctx);
	  if (print)
	    {
	      flint_printf("A(%wd) = ", delta);
	      fmpq_mpoly_print_pretty(test, (const char**) vars, ctx);
	      flint_printf("\n");
	    }
	  
	  humbert_get_mpoly(test, (const char**) vars, "B1", delta, ctx);
	  if (print)
	    {
	      flint_printf("B1(%wd) = ", delta);
	      fmpq_mpoly_print_pretty(test, (const char**) vars, ctx);
	      flint_printf("\n");
	    }
	  
	  humbert_get_mpoly(test, (const char**) vars, "B", delta, ctx);
	  if (print)
	    {
	      flint_printf("B(%wd) = ", delta);
	      fmpq_mpoly_print_pretty(test, (const char**) vars, ctx);
	      flint_printf("\n");
	    }
	  
	  humbert_get_mpoly(test, (const char**) vars, "B2", delta, ctx);
	  if (print)
	    {
	      flint_printf("B2(%wd) = ", delta);
	      fmpq_mpoly_print_pretty(test, (const char**) vars, ctx);
	      flint_printf("\n");
	    }
	}
    }

  humbert_vars_clear(vars);
  fmpq_mpoly_clear(test, ctx);
  fmpq_mpoly_ctx_clear(ctx);

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
