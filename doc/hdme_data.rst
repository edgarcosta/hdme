
**hdme_data.h** -- Accessing data stored in the hdme_data folder
================================================================

The **hdme_data** module is about accessing certain precomputed
multivariate polynomials stored as part of the library source
code. Each polynomial is stored in a one-line file.

The precomputed polynomials are stored in different folders:

* **hdme_data/igusa**: expressions of Igusa scalar covariants *I2, I4,
  I6', I10* in terms of curve coefficients; expressions of Mestre's
  scalar covariants *A,B,C,D* in terms of *I2, I4, I6', I10*.
  
* **hdme_data/mestre**: expressions of the coefficients of the cubic
  and conic in Mestre's algorithm in terms of the scalar covariants
  *A, B, C, D*.
  
* **hdme_data/humbert**: parametrizations of Humbert surfaces for all
  fundamental discriminants <=100, as computed in the following
  paper. N. D. Elkies, A. Kumar, K3 surfaces and equations for Hilbert
  modular surfaces, Algebra & Number Theory 8(10),
  pp. 2297-2441, 2014.

* **hdme_data/hilbert**: parametrizations of Hilbert surfaces for all
  fundamental discriminants <=20, found in the same paper.

* **hdme_data/gundlach**: conversion between Igusa scalar covariants
  and Gundlach invariants for discriminant 5.

  
Types and constants
-------------------

We make use of the Flint **fmpq_mpoly_t** type.

::
   HDME_DATA_STR_LEN

Upper bound on the length of each of these one-line files.

::
   HDME_DATA_FILE_LEN

Upper bound on the length of each file name.

::
   HDME_DATA_VAR_LEN
   

Upper bound on the length of each variable name.

::
   HDME_DATA_PATH

Constant string indicating the path to the **hdme_data** folder. This
uses **HDME_PATH**, a.k.a. **CURDIR** as defined in the Makefile.


Functions
---------

::
   char** hdme_data_vars_init(slong nb);

Return an initialized vector *vars* of *nb* strings of length
**HDME_DATA_VAR_LEN**.

::
   void hdme_data_vars_clear(char** vars, slong nb);

Free the memory used by *vars*.

::
   void hdme_data_vars_set(char** vars, const char* name, slong k);
   

Copy the string *name* to *vars[k]*.

::
   void hdme_data_read(fmpq_mpoly_t pol, const char** vars, const char* name,
		    const fmpq_mpoly_ctx_t ctx);
		    
Parse the content of the file named *name* as a multivariate
polynomial in the variables listed in *vars*, and store the result in
*pol*. The number of variables is encoded as part of the *ctx*
context. This function aborts on failure.

::
   void hdme_data_evaluate_acb(acb_t ev, const fmpq_mpoly_t pol, acb_srcptr vals,
			    const fmpq_mpoly_ctx_t ctx, slong prec);

Set *ev* to the evaluation of *pol* after substitution of all
variables by the values stored in *vals*. The number of variables is
encoded as part of the *ctx* context.
