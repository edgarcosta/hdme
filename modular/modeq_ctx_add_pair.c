
#include "modular.h"

void modeq_ctx_add_pair(modeq_ctx_t ctx, slong i1, slong i2)
{
  slong k = modeq_ctx_nb_pairs(ctx);
  
  if (2*(k+1) > modeq_ctx_alloc_pairs(ctx))
    {
      modeq_ctx_alloc_pairs(ctx) += MODEQ_CTX_ALLOC;
      modeq_ctx_pairs(ctx) = flint_realloc(modeq_ctx_pairs(ctx),
					   modeq_ctx_alloc_pairs(ctx));
    }
  
  modeq_ctx_nb_pairs(ctx) = k+1;
  modeq_ctx_pair(ctx, k)[0] = i1;
  modeq_ctx_pair(ctx, k)[1] = i2;    
}
