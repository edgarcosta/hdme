
#include "modular.h"

void siegel_modeq_rescale(fmpz_t scal, fmpq* j, slong ell)
{
  fmpz_t den;
  fmpz* num;
  fmpq* j_mod;
  slong exp = (5 * siegel_nb_cosets(ell))/3;

  fmpz_init(den);
  num = _fmpz_vec_init(3);
  j_mod = _fmpq_vec_init(3);

  fmpq_mul(&j_mod[0], &j[0], &j[0]);
  fmpq_set(&j_mod[1], &j[1]);
  fmpq_set(&j_mod[2], &j[2]);
  modeq_input_get_fmpz(den, num, j_mod, 3);
  fmpz_pow_ui(scal, den, exp);

  fmpz_clear(den);
  _fmpz_vec_clear(num, 3);
  _fmpq_vec_clear(j_mod, 3);
}
