
#include "siegel.h"

static void
siegel_test_matrix_g1(sp2gz_t u, slong j)
{
  switch(j)
    {
    case 0:
      sp2gz_J(u);
      break;
    default:
      flint_printf("Invalid index in siegel_test_matrix: %d\n", j);
      flint_abort();
    }
}

static void
siegel_test_matrix_g2(sp2gz_t u, slong j)
{
  fmpz_mat_zero(&u->a);
  fmpz_mat_zero(&u->b);
  fmpz_mat_zero(&u->c);
  fmpz_mat_zero(&u->d);

  if (j < 15)
    {
      fmpz_mat_zero(&u->a);
      fmpz_mat_one(&u->c);
      fmpz_mat_neg(&u->b, &u->c);
    }
  if (15 <= j && j < 17)
    {
      fmpz_mat_one(&u->a);
      fmpz_mat_neg(&u->b, &u->a);
    }
  if (17 <= j && j < 19)
    {
      fmpz_mat_zero(&u->b);
      fmpz_one(fmpz_mat_entry(&u->c, 0, 0));
      fmpz_set_si(fmpz_mat_entry(&u->c, 0, 1), -1);
      fmpz_set_si(fmpz_mat_entry(&u->c, 1, 0), -1);
      fmpz_one(fmpz_mat_entry(&u->c, 1, 1));
    }
  
  switch(j)
    {
    case 0:
      fmpz_mat_zero(&u->d);
      break;
    case 1:
      fmpz_one(fmpz_mat_entry(&u->d, 0, 0));
      break;
    case 2:
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 0), -1);
      break;
    case 3:
      fmpz_one(fmpz_mat_entry(&u->d, 1, 1));
      break;
    case 4:
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 1), -1);
      break;
    case 5:
      fmpz_mat_one(&u->d);
      break;
    case 6:
      fmpz_mat_one(&u->d);
      fmpz_mat_neg(&u->d, &u->d);
      break;
    case 7:
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 0), -1);
      fmpz_one(fmpz_mat_entry(&u->d, 1, 1));
      break;
    case 8:
      fmpz_one(fmpz_mat_entry(&u->d, 0, 0));
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 1), -1);
      break;
    case 9:
      fmpz_one(fmpz_mat_entry(&u->d, 0, 1));
      fmpz_one(fmpz_mat_entry(&u->d, 1, 0));
      break;
    case 10:
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 1), -1);
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 0), -1);
      break;
    case 11:
      fmpz_one(fmpz_mat_entry(&u->d, 0, 0));
      fmpz_one(fmpz_mat_entry(&u->d, 0, 1));
      fmpz_one(fmpz_mat_entry(&u->d, 1, 0));
      break;
    case 12:
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 0), -1);
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 1), -1);
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 0), -1);
      break;
    case 13:
      fmpz_one(fmpz_mat_entry(&u->d, 0, 1));
      fmpz_one(fmpz_mat_entry(&u->d, 1, 0));
      fmpz_one(fmpz_mat_entry(&u->d, 1, 1));
      break;
    case 14:
      fmpz_set_si(fmpz_mat_entry(&u->d, 0, 1), -1);
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 0), -1);
      fmpz_set_si(fmpz_mat_entry(&u->d, 1, 1), -1);
      break;
    case 15:
      fmpz_one(fmpz_mat_entry(&u->c, 0, 0));
      fmpz_one(fmpz_mat_entry(&u->d, 1, 1));
      break;
    case 16:
      fmpz_one(fmpz_mat_entry(&u->c, 1, 1));
      fmpz_one(fmpz_mat_entry(&u->d, 0, 0));
      break;
    case 17:
      fmpz_mat_one(&u->a);
      fmpz_mat_one(&u->d);
      break;
    case 18:
      fmpz_mat_one(&u->a);
      fmpz_mat_neg(&u->a, &u->a);
      fmpz_mat_set(&u->d, &u->a);
      break;
    default:
      flint_printf("Invalid index in siegel_test_matrix: %d\n", j);
      flint_abort();
    }
}

void
siegel_test_matrix(sp2gz_t u, slong j)
{
  slong g = u->g;
  switch(g)
    {
    case 1:
      siegel_test_matrix_g1(u, j);
      break;
    case 2:
      siegel_test_matrix_g2(u, j);
      break;
    default:
      flint_printf("siegel_test_matrix not implemented for g = %d\n", g);
      flint_abort();
    }
}
