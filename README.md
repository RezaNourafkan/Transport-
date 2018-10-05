# Transport-
Longitudinal transport polarization
The trasport code is a Fortran code for calculating the longitudinal transport polarization, assuming a momentum independent self-energy (DMFT approximation). More especifiaclly, it calculates the paramagnetic contribution of Eq. (11) at the following paper

arXiv:1807.03852

Analytic continuation is required to obtain optical conductivity. A python code for Pade analytic continuation is accompanied. 

Instalation:
Go to src directory, revise the makefile by providing the path to lapack and blas library  and command
make

It will create an exectable file, called "Tran", outside of src file.

Input:


