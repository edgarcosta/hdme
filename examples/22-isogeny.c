
#include "modular.h"

int main()
{
  fmpq* j;
  fmpq* all_isog_j;

  slong ell = 2;
  slong max_nb_roots = siegel_nb_cosets(ell);
  slong nb_roots = 0;
  slong i;

  j = _fmpq_vec_init(3);
  all_isog_j = _fmpq_vec_init(3*max_nb_roots);

  fmpq_set_str(&j[0], "-140575574447/84500", 10);
  fmpq_set_str(&j[1], "1571037533912/21125", 10);
  fmpq_set_str(&j[2], "31651191767790319272856/446265625", 10);

  siegel_modeq_isog_invariants_Q(&nb_roots, all_isog_j, j, ell);
  for (i = 0; i < 3*nb_roots; i++)
    {
      fmpq_print(&all_isog_j[i]); flint_printf("\n");
    }
    
  
  _fmpq_vec_clear(j, 3);
  _fmpq_vec_clear(all_isog_j, 3*max_nb_roots);
  return EXIT_SUCCESS;
}
