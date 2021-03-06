.. Wright documentation master file, created by
   sphinx-quickstart on Fri Nov 26 08:12:18 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Wright's documentation!
==================================

This small library contains implementation in C, Fortran and MATLAB for computing
the Wright function of the second kind, often called Mainardi's function, on the
real line.

How to obtain
-------------

The latest version of this code can be obtained from GitHub

.. code-block:: none

   git clone git@github.com:Cirdans-Home/mwright.git

How to compile
--------------
The code can be compiled using the CMake toolchain in the usual way.

.. code-block:: none

   mkdir build
   cd build
   ccmake ..
   make

To build the :code:`c` routines you'll need a working version of the
`ARB library <https://fredrikj.net/arb/>`_, this is a a C library for
arbitrary-precision floating-point ball arithmetic. It is used here to
implement the evaluation of the Wright function through its series
expansion. We use it to measure the error with respect to the algorithm
based on the inversion of the Laplace transform.

.. warning::
   If you want to build this documentation you'll need also a working version
   of Sphinx, together with the extensions :code:`'sphinx.ext.autodoc'`,
   :code:`'sphinx_c_autodoc'`, :code:`'sphinx.ext.mathjax'`,
   :code:`'sphinxfortran.fortran_domain'`, :code:`'sphinxfortran.fortran_autodoc'`,
   and :code:`sphinxcontrib.matlab`.
   This step **can be omitted** from the :code:`ccmake` settings.

How to cite
-----------

Please, if you use this code for a scientific publication cite the work:

- L. Aceto and F. Durastante. “Efficient computation of the Wright function and its applications to fractional diffusion-wave equations”. 2022. arXiv: `2202.00397 <https://arxiv.org/abs/2202.00397>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   ccode
   fcode
   mcode
   testprograms



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
