
#include "modular.h"

void hilbert_modeq_gundlach_rescale(fmpz_t scal, fmpq* g, slong ell, slong delta)
{
  fmpz_t den;
  fmpz* num;
  fmpq* g_mod;
  slong wt = 2*10*hilbert_nb_cosets(ell, delta);

  fmpz_init(den);
  num = _fmpz_vec_init(2);
  g_mod = _fmpq_vec_init(2);

  fmpq_pow_si(&g_mod[0], &g[0], wt/6);
  fmpq_pow_si(&g_mod[1], &g[1], wt/6);
  modeq_input_get_fmpz(den, num, g_mod, 2);
  fmpz_pow_ui(scal, den, 1);

  fmpz_clear(den);
  _fmpz_vec_clear(num, 2);
  _fmpq_vec_clear(g_mod, 2);
}
