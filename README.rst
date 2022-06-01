Description
===========

This repository contains a C library for evaluating higher-dimensional
modular equations by analytic methods. As it stands now, only a few
types of modular equations in genus g=2 are supported, either over the
rationals or prime finite fields:

- Modular equations of Siegel type in Igusa invariants;

- Modular equations of Hilbert type in Igusa invariants, using
  parametrizations of Humbert surfaces for all fundamental
  discriminants less than 100;

- Modular equations of Hilbert type in Gundlach invariants for
  discriminant 5.

Along the way, the library provides functionality for the evaluation
of genus 2 theta constants in quasi-linear time, and for the action of
Hecke operators on Siegel modular forms. The algorithms are described
in the associated paper `[K21]`_.

Prerequisites
=============

An installation of `Arb`_ version 2.21.0 or later, and `Flint`_
version 2.8.0 or later.

Installation
============

Run
::
  ./configure
  
to create the Makefile. To compile the library, run
::
  make
  
This will produce a shared library in the source directory. The
Makefile also responds to the standard commands ::
  make tests
  make check
  make clean
  
Documentation
=============

A detailed documentation for each function existed in previous
versions, but is now outdated. For examples of usage, look at the test
files in the `test` subfolder of each module.

Credits
=======

- The code for fast evaluation of theta constants using Arb is
  inspired by the `cmh`_ library by Andreas Enge and Emmanuel Thomé.

- The code written by `Enea Milio`_ has been useful in the
  implementation of Mestre's algorithm and to implement the inversion
  of the Hilbert embedding.

- The design and semantics of many functions are inspired from existing
  functions in the libraries Flint and Arb, maintained by William Hart
  and Fredrik Johansson.
  
.. _[K21]: https://arxiv.org/abs/2010.10094
.. _Flint: https://flintlib.org
.. _Arb: https://arblib.org
.. _cmh: https://gitlab.inria.fr/cmh/cmh
.. _Enea Milio: https://members.loria.fr/EMilio/modular-polynomials/
