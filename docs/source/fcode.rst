Fortran Routines
================

This is the documentation of the Fortran routines for calculating the Wright
function. Unlike their counterpart in C, they use ordinary arithmetic with
different types of precision (32,64,128) for the calculation.

We have implemented both the scalar version of the function, and the vectorized
counterparts with respect to the :math:`x` variable. They are all
*private* inside the module, and are accesible from a unique interface called
`wright`.

.. f:autosrcfile:: wright.f90


.. toctree::
   :maxdepth: 2
   :caption: Contents:
