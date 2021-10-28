Description
===========

This repository contains a C library for evaluating higher-dimensional
modular equations by analytic methods. As it stands now, only a few
types of modular equations in genus g=2 are supported:

- Modular equations of Siegel type in Igusa invariants;

- Modular equations of Hilbert type in Igusa invariants, using
  parametrizations of Humbert surfaces for all fundamental
  discriminants less than 100;

- Modular equations of Hilbert type in Gundlach invariants for
  discriminant 5.

Along the way, the library provides functionality for the evaluation
of genus 2 theta constants in quasi-linear time. The algorithms are
described in the associated paper `[K21]`_.

Prerequisites
=============

An installation of `Arb`_ version 2.21.0 or later, and `Flint`_
version 2.8.0 or later.

Installation
============

Save the file `Makefile.in` as `Makefile` at the top of the source
tree, and modify the variables `CURDIR` and `MYDIR` as appropriate. Then, run::
  make
  
to compile the library. This will produce an archive file named::
  libhdme.a
  
in the source directory. The Makefile also responds to the standard commands::
  make tests
  make check
  make clean

Documentation
=============

See the `doc` folder for the documentation of the individual modules.

Credits
=======

The fast evaluation of theta constants in Arb is inspired by the
`cmh`_ library by Andreas Enge and Emmanuel Thomé.

The code written by `Enea Milio`_ has been useful in the
implementation of Mestre's algorithm.

The various parametrizations of Hilbert and Humbert surfaces used here
appear in the paper `[EK14]`_ by Noam D. Elkies and Abhinav Kumar.
