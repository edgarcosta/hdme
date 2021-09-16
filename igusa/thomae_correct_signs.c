
#include "igusa.h"

int thomae_correct_signs(slong* perm, slong* signs, acb_srcptr roots,
			 acb_srcptr I, slong prec)
{
  slong current_prec = thomae_startprec(prec);
  slong nb_candidates = 720 * 16;
  slong p; /* 0 to 719 */
  slong s; /* 0 to 15 */
  slong correct_perm = -1;
  slong correct_signs = -1;
  acb_ptr new_roots;
  acb_ptr ros;
  acb_ptr th4, th2;
  acb_ptr j, j_test;
  slong k;
  acb_mat_t tau;
  int res = 1;
  int tau_success;
  int v = THOMAE_VERBOSE;
  
  int removed_ps[720][16];
  int removed_p[720];

  new_roots = _acb_vec_init(6);
  ros = _acb_vec_init(3);
  th4 = _acb_vec_init(16);
  th2 = _acb_vec_init(16);
  j = _acb_vec_init(3);
  j_test = _acb_vec_init(3);
  acb_mat_init(tau, 2, 2);

  for (p = 0; p < 720; p++)
    {
      removed_p[p] = 0;
      for (s = 0; s < 16; s++) removed_ps[p][s] = 0;
    }
  
  while ((nb_candidates > 1) && (current_prec < prec))
    {
      if(v) flint_printf("(thomae_correct_signs) Trying precision %wd\n", current_prec);
      nb_candidates = 0;
      for (p = 0; p < 720; p++)
	{
	  if (!removed_p[p])
	    {
	      removed_p[p] = 1;
	      thomae_reorder(new_roots, roots, p);
	      thomae_rosenhain(ros, new_roots, current_prec);
	      thomae_theta4(th4, ros, current_prec);
	      for (s = 0; s < 16; s++)
		{
		  if (!removed_ps[p][s])
		    {
		      removed_ps[p][s] = 1;
		      thomae_theta2(th2, th4, ros, s, current_prec);
		      if (!thomae_discard(th2, current_prec))
			{
			  /* flint_printf("Keeping (p,s) = (%wd,%wd)\n", p, s);
			     for (k = 0; k < 4; k++)
			     {
			     acb_printd(&th2[k], 10); flint_printf("\n");
			     } */
			  tau_success = theta2_inverse(tau, th2, current_prec);
			  /* We have to keep them in case of failure to compute. */
			  if (!tau_success || thomae_keep_candidate(tau, I, current_prec))
			    {
			      correct_perm = p;
			      correct_signs = s;
			      removed_p[p] = 0;
			      removed_ps[p][s] = 0;
			      nb_candidates++;
			      /* acb_mat_printd(tau, 10); 
				 flint_printf("Succeeded computing tau: %d\n", tau_success); */
			    }
			  /* Else, tau was discarded: do nothing */
			}
		      /* Else, th2 was discarded: do nothing */
		    }
		  /* Else (s,p) was removed: do nothing */
		}
	      /* End for(s) loop */
	    }
	  /* Else p was removed: do nothing */
	}
      /* End for(p) loop */
      if(v) flint_printf("(thomae_correct_signs) Remaining candidates: %wd\n", nb_candidates);
      current_prec *= THOMAE_MULPREC;
    }
  
  /* Then we have a last run at prec: this time we want to succeed in computing tau */
  if (nb_candidates > 1)
    {
      if(v) flint_printf("(thomae_correct_signs) Last run at precision %wd\n", prec);
      nb_candidates = 0;
      igusa_from_cov(j_test, I, prec);
      for (p = 0; p < 720; p++)
	{
	  if (!removed_p[p])
	    {
	      removed_p[p] = 1;
	      thomae_reorder(new_roots, roots, p);
	      thomae_rosenhain(ros, new_roots, prec);
	      thomae_theta4(th4, ros, prec);
	      for (s = 0; s < 16; s++)
		{
		  if (!removed_ps[p][s])
		    {
		      removed_ps[p][s] = 1;
		      thomae_theta2(th2, th4, ros, s, current_prec);
		      if (!thomae_discard(th2, current_prec))
			{
			  tau_success = theta2_inverse(tau, th2, current_prec);
			  if (tau_success) tau_success = igusa_from_tau(j, tau, prec);
			  if (tau_success)
			    {
			      for (k = 0; k < 3; k++)
				{
				  if (!acb_overlaps(&j[k], &j_test[k])) tau_success = 0;
				}
			    }
			  if (tau_success)
			    {
			      correct_perm = p;
			      correct_signs = s;
			      nb_candidates++;
			      
			      *perm = correct_perm;
			      *signs = correct_signs;
			      /* acb_mat_printd(tau, 10); */
			      goto exit;
			    }
			}
		    }
		}
	    }
	}
      if(v) flint_printf("(thomae_correct_signs) Remaining candidates: %wd\n", nb_candidates);
    }
  if (nb_candidates == 0) res = 0;
  else
    {
      *perm = correct_perm;
      *signs = correct_signs;
    }
  goto exit;

 exit:
  {
    _acb_vec_clear(new_roots, 6);
    _acb_vec_clear(ros, 3);
    _acb_vec_clear(th4, 16);
    _acb_vec_clear(th2, 16);
    _acb_vec_clear(j, 3);
    _acb_vec_clear(j_test, 3);
    acb_mat_clear(tau);
    return res;
  }
}
