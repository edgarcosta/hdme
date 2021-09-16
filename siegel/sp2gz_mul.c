
#include "siegel.h"

void
sp2gz_mul(sp2gz_t r, const sp2gz_t m, const sp2gz_t n)
{
  fmpz_mat_t block_product;
  sp2gz_t t;
  
  if (r == m || r == n)
    {
      sp2gz_init(t, m->g);
      sp2gz_mul(t, m, n);
      sp2gz_swap(t, r);
      sp2gz_clear(t);
      return;
    }
  
  fmpz_mat_init(block_product, m->g, m->g);
  
  fmpz_mat_mul(&r->a, &m->a, &n->a);
  fmpz_mat_mul(block_product, &m->b, &n->c);
  fmpz_mat_add(&r->a, &r->a, block_product);
  
  fmpz_mat_mul(&r->b, &m->a, &n->b);
  fmpz_mat_mul(block_product, &m->b, &n->d);
  fmpz_mat_add(&r->b, &r->b, block_product);
  
  fmpz_mat_mul(&r->c, &m->c, &n->a);
  fmpz_mat_mul(block_product, &m->d, &n->c);
  fmpz_mat_add(&r->c, &r->c, block_product);
  
  fmpz_mat_mul(&r->d, &m->c, &n->b);
  fmpz_mat_mul(block_product, &m->d, &n->d);
  fmpz_mat_add(&r->d, &r->d, block_product);
  
  fmpz_mat_clear(block_product);
}
    
