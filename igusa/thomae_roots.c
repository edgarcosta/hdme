
#include "igusa.h"

/* Use the generic acb root finding algorithm; then check all roots
   are well isolated */

/* Max number of iterations of the numeric method: 0 does not work,
   use prec */

int thomae_roots(acb_ptr roots, const acb_poly_t crv, slong prec)
{
  slong n;
  int res;
  n = acb_poly_find_roots(roots, crv, NULL, prec, prec);
  res = (n == acb_poly_degree(crv));
  if (!res) {
    if (get_thomae_verbose()) {
      flint_printf("(thomae_roots) Warning: unable to compute roots\n");
      flint_printf("(thomae_roots) Curve equation:\n");
      acb_poly_printd(crv, 10);
      flint_printf("\n(thomae_roots) Number of isolated roots: %wd\n", n);
      flint_printf("(thomae_roots) Working precision: %wd\n", prec);
    }
  }
  return res;
}
