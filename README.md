# Transport-
**Longitudinal transport polarization**

The trasport code is a Fortran code for calculating the longitudinal transport polarization, assuming a momentum independent self-energy (DMFT approximation). More especifiaclly, it calculates the paramagnetic contribution of Eq. (11) at the following paper

**arXiv:1807.03852**

Analytic continuation is required to obtain optical conductivity. A python (2.7) code (called Pade.py) for Pade analytic continuation is accompanied. It can be used as follows:

python Pade.py  OpticalPolarization.dat -n 2  -nfreq  64  -nene 1600 -enecutoff 12 -o PadeOut.dat

where

-n is column number for data in OpticalPolarization.dat file

-nfreq is number of bosonic Matsubara frequency used in analytic continuation

-nene is number of energy grid on real axis

-enecutoff is energy cut off

-o is the output file name

**The optical conductivity is given by third column of PadeOut.dat divided by first column.**

Instalation:
Go to src directory, revise the makefile by providing the path to lapack and blas library  and command
make

It will create an exectable file, called "Tran", outside of src file.

**Input:**
The iput file should be saved in a directory called "Input" and includes:

Input.dat : a text file including some parameters

case.kpoints : a text file including the reciprocal lattice vectors and number of k points

Case.Hk : A binary file including Hamiltonian in momementum space. It should be written as follows in Fortran:

      write(NUnitOut) Nnd            ! Number of Orbitals
      write(NUnitOut) Numkpoints     ! Number of K points
      do kindex = 1, Numkpoints
        write(NUnitOut) ((Hk(norb1,norb2), norb2=1,Nnd), norb1=1,Nnd) ! Hk is Hamiltonian matrix
      end do
      
Case.vk : A binary file including current vertex in momementum space. It should be written as follows:

      do kindex=1, Numkpoints
        write(NUnitOut) kindex, Nnd
        do norb1 = 1, Nnd
          do norb2 = 1, Nnd
            write(NUnitOut) norb1,norb2, (Mat(n,norb1,norb2,kindex), n=1, 3)
          end do
        end do  
        
 Case.SelfEnIA : Momentum-independent Self energy on imaginary axis written as
 
     do N = NiwnMin, NiwnMax-1
       Omega_n = (Two*n + One)*Pi_d/InverseTemperature
       write(NUnitOut,formmat,advance='no') N, Omega_n, &
         (((SelfEnergy(Ns,Norb1,Norb2,N), NOrb2=1,Nnd),NOrb1=1,Nnd), ns=1,Nspin)
     end do
 
 
 **Outputs:** 
 
 Local Green's function
 
 Optical polarization
 
 


