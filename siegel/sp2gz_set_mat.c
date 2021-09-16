
#include "siegel.h"

void
sp2gz_set_mat(sp2gz_t m, const fmpz_mat_t n)
{
  fmpz_mat_t window;
  slong g = m->g;

  fmpz_mat_window_init(window, n, 0, 0, g, g);
  fmpz_mat_set(&m->a, window);
  fmpz_mat_window_clear(window);

  fmpz_mat_window_init(window, n, 0, g, g, 2*g);
  fmpz_mat_set(&m->b, window);
  fmpz_mat_window_clear(window);

  fmpz_mat_window_init(window, n, g, 0, 2*g, g);
  fmpz_mat_set(&m->c, window);
  fmpz_mat_window_clear(window);

  fmpz_mat_window_init(window, n, g, g, 2*g, 2*g);
  fmpz_mat_set(&m->d, window);
  fmpz_mat_window_clear(window);
}
