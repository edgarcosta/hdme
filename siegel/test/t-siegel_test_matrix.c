
#include "siegel.h"

int main()
{
  slong i, j, g;
  sp2gz_t ui, uj;
  
  flint_printf("siegel_test_matrix....");
  fflush(stdout);

  /* Test matrices are distinct and symplectic */
  for (g = 1; g <= 2; g++)
    {
      sp2gz_init(ui, g);
      sp2gz_init(uj, g);
      for (i = 0; i < siegel_nb_test_matrices(g); i++)
	{
	  siegel_test_matrix(ui, i);
	  /* sp2gz_print(ui); flint_printf("\n\n"); */
	  if (!sp2gz_is_correct(ui))
	    {
	      flint_printf("FAIL (not symplectic)\n");
	      flint_printf("i = %wd\n", i);
	      sp2gz_print(ui); flint_printf("\n\n");
	      flint_abort();
	    }
	  for (j = 0; j < siegel_nb_test_matrices(g); j++)
	    {
	      siegel_test_matrix(uj, j);
	      if (i != j && sp2gz_equal(ui, uj))
		{
		  flint_printf("FAIL (equal matrices)\n");
		  flint_printf("(i, j) = (%wd, %wd)\n", i, j);
		  sp2gz_print(ui); flint_printf("\n\n");
		  flint_abort();
		}
	    }
	}
      sp2gz_clear(ui);
      sp2gz_clear(uj);
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
