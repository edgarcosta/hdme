
#include "igusa.h"

void cardona(acb_poly_t crv, acb_srcptr IC, slong prec)
{
  acb_ptr ABCD;
  acb_ptr Aij; /* 11, 12, 22, 33 */
  acb_ptr aijk; /* 111, 112, 122, 133, 222, 233 */
  acb_poly_t P1, P2, P3;
  acb_t c;
  acb_poly_t term;
  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;

  ABCD = _acb_vec_init(4);
  Aij = _acb_vec_init(4);
  aijk = _acb_vec_init(6);
  acb_poly_init(P1);
  acb_poly_init(P2);
  acb_poly_init(P3);
  acb_init(c);
  acb_poly_init(term);
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  hdme_data_vars_set(vars, "A", 0);
  hdme_data_vars_set(vars, "B", 1);
  hdme_data_vars_set(vars, "C", 2);
  hdme_data_vars_set(vars, "D", 3);

  /* Set coefficients */
  igusa_ABCD_from_IC(ABCD, IC, prec);
  
  hdme_data_read(pol, (const char**) vars, "cardona/A11", ctx);
  hdme_data_evaluate_acb(&Aij[0], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A12", ctx);
  hdme_data_evaluate_acb(&Aij[1], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A22", ctx);
  hdme_data_evaluate_acb(&Aij[2], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/A33", ctx);
  hdme_data_evaluate_acb(&Aij[3], pol, ABCD, ctx, prec);

  hdme_data_read(pol, (const char**) vars, "cardona/a111", ctx);
  hdme_data_evaluate_acb(&aijk[0], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a112", ctx);
  hdme_data_evaluate_acb(&aijk[1], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a122", ctx);
  hdme_data_evaluate_acb(&aijk[2], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a133", ctx);
  hdme_data_evaluate_acb(&aijk[3], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a222", ctx);
  hdme_data_evaluate_acb(&aijk[4], pol, ABCD, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "cardona/a233", ctx);
  hdme_data_evaluate_acb(&aijk[5], pol, ABCD, ctx, prec);

  /* Set P1, P2, P3 */
  acb_mul_si(c, &Aij[1], -2, prec);
  acb_poly_set_coeff_acb(P1, 0, c);
  acb_mul_si(c, &Aij[2], -2, prec);
  acb_poly_set_coeff_acb(P2, 1, c);

  acb_poly_set_coeff_acb(P2, 0, &Aij[0]);
  acb_neg(c, &Aij[2]);
  acb_poly_set_coeff_acb(P2, 2, c);

  acb_poly_set_coeff_acb(P3, 0, &Aij[0]);
  acb_mul_si(c, &Aij[1], 2, prec);
  acb_poly_set_coeff_acb(P3, 1, c);
  acb_poly_set_coeff_acb(P3, 2, &Aij[2]);

  /* Set crv */
  acb_poly_zero(crv);
  
  acb_poly_pow_ui(term, P1, 3, prec);
  acb_mul(c, &Aij[3], &aijk[0], prec);
  acb_neg(c, c);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);

  acb_poly_pow_ui(term, P1, 2, prec);
  acb_poly_mul(term, term, P2, prec);
  acb_mul(c, &Aij[3], &aijk[1], prec);
  acb_mul_si(c, c, -3, prec);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);

  acb_poly_pow_ui(term, P2, 2, prec);
  acb_poly_mul(term, term, P1, prec);
  acb_mul(c, &Aij[3], &aijk[2], prec);
  acb_mul_si(c, c, -3, prec);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P3, 2, prec);
  acb_poly_mul(term, term, P1, prec);
  acb_mul(c, &Aij[2], &aijk[3], prec);
  acb_mul_si(c, c, 3, prec);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P3, 3, prec);
  acb_mul(c, &Aij[3], &aijk[4], prec);
  acb_neg(c, c);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);
  
  acb_poly_pow_ui(term, P3, 2, prec);
  acb_poly_mul(term, term, P2, prec);
  acb_mul(c, &Aij[2], &aijk[5], prec);
  acb_mul_si(c, c, 3, prec);
  acb_poly_scalar_mul_acb(term, term, c, prec);
  acb_poly_add(crv, crv, term, prec);

  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  
  _acb_vec_clear(ABCD, 4);
  _acb_vec_clear(Aij, 4);
  _acb_vec_clear(aijk, 6);
  acb_poly_clear(P1);
  acb_poly_clear(P2);
  acb_poly_clear(P3);
  acb_clear(c);
  acb_poly_clear(term);
}
