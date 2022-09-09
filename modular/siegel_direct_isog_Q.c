
#include "modular.h"

int siegel_direct_isog_Q(slong* nb, fmpz* all_I, fmpz* I, slong ell)
{
  slong prec = siegel_modeq_startprec(I, ell);
  hecke_t H;
  int stop = 0;
  int res = 1;
  *nb = 0; /* we might never call hecke_all_isog_Q */

  hecke_init(H, siegel_nb_cosets(ell));

  while (!stop) {
    res = hecke_set_I_fmpz(H, I, prec);

    if (res) res = hecke_collect_siegel(H, ell, prec);
    if (res) hecke_make_integral(H, I, prec);

    if (res && hecke_has_integral_precision(H, prec)) {
      res = hecke_all_isog_Q(nb, all_I, H, I, prec);
      if (res) stop = 1;
    }

    prec = modeq_nextprec_generic(prec);
    stop = modeq_stop(res, prec);
  }

  hecke_clear(H);
  return res;
}
