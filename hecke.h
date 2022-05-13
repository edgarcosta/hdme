/*
    Copyright (C) 2022 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#ifndef HECKE_H
#define HECKE_H

#include <acb.h>
#include <acb_mat.h>
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

#include "hilbert.h"
#include "igusa.h"
#include "siegel.h"
#include "theta.h"

#define HECKE_RED_TOL_BITS 50

#ifndef HECKE_VERBOSE
#define HECKE_VERBOSE 0
#endif

/* Define a data structure hecke_t that contains all the necessary
   data related to Hecke correspondences over the complex numbers,
   i.e. data about all period matrices isogenous to a given one. */

typedef struct
{
  /* Context information */
  acb_mat_t tau; /* Original period matrix */
  slong ell; /* Type of isogenies */
  acb_ptr t1t2; /* Original periods (Hilbert case only) */
  fmpz_poly_t beta; /* Type of isogenies (Hilbert case only) */

  /* Data to be computed */
  slong nb; /* Number of isogenous matrices */
  fmpz_mat_struct* cosets; /* General symplectic matrices acting */
  acb_mat_struct* isog; /* Isogenous period matrices in fundamental domain */
  acb_mat_struct* stars; /* Cocycles C*tau+D */
  acb_ptr stardets; /* Determinants of stars */
  acb_ptr theta2; /* Values of theta constants at isog */
  acb_ptr I; /* Values of Igusa-Clebsch covariants at isog */    
} hecke_struct;

typedef hecke_struct hecke_t[1];

#define hecke_tau(H) ((H)->tau)
#define hecke_ell(H) ((H)->ell)
#define hecke_t1t2(H) ((H)->t1t2)
#define hecke_beta(H) ((H)->beta)
#define hecke_nb(H) ((H)->nb)

#define hecke_coset(H, k) (&(H)->cosets[(k)])
#define hecke_isog(H, k) (&(H)->isog[(k)])
#define hecke_star(H, k) (&(H)->stars[(k)])
#define hecke_stardet(H, k) (&(H)->stardets[(k)])
#define hecke_theta2(H, k) (&(H)->theta2[16*(k)])
#define hecke_I(H, k) (&(H)->I[4*(k)])

void hecke_init(hecke_t H, slong nb);

void hecke_clear(hecke_t H);

int hecke_set_entry(hecke_t H, slong k, const fmpz_mat_t gamma, slong prec);


slong siegel_nb_cosets(slong ell);

void siegel_coset(fmpz_mat_t m, slong k, slong ell);

int hecke_set_siegel(hecke_t H, const acb_mat_t tau, slong ell, slong prec);


slong hilbert_nb_cosets(slong ell, slong delta);

void hilbert_coset(fmpz_poly_mat_t m, slong k, slong ell, slong delta);

int hecke_set_hilbert(hecke_t H, acb_srcptr t, const fmpz_poly_t beta,
		      slong ell, slong delta, slong prec);

int hecke_set_hilbert_sym(hecke_t H, acb_srcptr t, const fmpz_poly_t beta,
			  slong ell, slong delta, slong prec);

#endif 

