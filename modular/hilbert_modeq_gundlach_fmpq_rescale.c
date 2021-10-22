
#include "modular.h"

void hilbert_modeq_gundlach_fmpq_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta)
{
  fmpz_t den;
  fmpz* num;
  fmpq* g_mod;
  slong exp = (2*hilbert_nb_cosets(ell))/3;

  fmpz_init(den);
  num = _fmpz_vec_init(2);
  g_mod = _fmpq_vec_init(2);

  fmpq_set(&g_mod[0], &g[0]);
  fmpq_pow_si(&g_mod[1], &g[1], 5);
  siegel_modeq_fmpz_input(den, num, g_mod, 2);
  fmpz_pow_ui(scal, den, exp);

  fmpz_clear(den);
  _fmpz_vec_clear(num, 2);
  _fmpq_vec_clear(j_mod, 2);
}
