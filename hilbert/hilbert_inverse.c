
#include "hilbert.h"

int hilbert_inverse(acb_t t1, acb_t t2, sp2gz_t eta, const acb_mat_t tau,
		    slong delta, slong prec)
{
  fmpz* abcde;

  abcde = _fmpz_vec_init(5);

  _fmpz_vec_clear(abcde, 5);
}
