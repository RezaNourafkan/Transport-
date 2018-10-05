#
#
# In the Pade approximation a function f(z)  in the complex plane, z, is described as the ratio between to polynomials P(z)  and Q(z),
# f(z)  = P(z)/Q(z) . The function is fitted to the output of a calculation on imaginary axis so that the results for certain imaginary 
# frequencies  N are reproduced exactly. The analytical continuation is then performed by evaluating the function on the real axis.

#!/usr/bin/env python
from types import *
import numpy
import string

def ReadDataFile(FileName, ColumnNumber):
    """ Reading data file """
    import os.path
    if not (os.path.exists(FileName)): raise IOError, "File %s does not exists"%FileName   # Check if file  exist
    FileHandle = open(FileName, "r")                     # open File
#   Blank lines and lines started with a number sign (#) in the data file will be ignored
    DataIn = numpy.loadtxt(FileHandle, unpack=True, usecols=[0,ColumnNumber-1, ColumnNumber])
    FileHandle.close()
    return numpy.transpose(DataIn)

def WriteDataFile(FileName, Data):
    """ Write data file """
    (NumRow, NumColumn) = numpy.shape(Data)
    format = " %14.10f " * NumColumn
    numpy.savetxt(FileName, Data, fmt=(format))
    return None

def Pade(EnergyMesh, Omega_m, GFniwn):
    """ Analytic continuation of a fermionic Greens function.
        The asymptotic behavior at large frequency MUST be like 1/(i\omega_m)    """
#
#   expanding the function in a continued fraction form, C_M(z), (this part is only works for fermionic GFns) 
#            c_1   c_2(z-z_1)       c_M(z-z_{M-1})
#   C_M(z) = ----  ----------- .... ---------------
#            1+        1+           1+(z-z_{M})
#   and compute the coefficients c_i.
#   M is number of Matsubara frequencies (so the degree of continued fraction is determined by number of frequencies)
#   Coefficients are determined in a way that C_M(z_i) = u_i, i=1, .., M, where u_i is the value of
#   the function in z_i. This is fulfilled if c_i satisfy a recursion relation
#   c_i = g_i(z_i), g_1(z_i) = u_i, i=1,2,3,...,M
#   g_i(z_j) = (g_{i-1}(z_{i-1})-g_{i-1}(z_j))/((z_j-z_{i-1})*g_{i-1}(z_j)),   i >= 2  
#
    NumFreq =  len(Omega_m)
    ContinuedFractionCoeff = numpy.zeros((NumFreq), dtype=numpy.complex)  
    ContinuedFractionCoeff[0] = GFniwn[0]
    for n in range(1,NumFreq):
        for m in range(n,NumFreq):
            GFniwn[m] = (GFniwn[n-1]-GFniwn[m])/(1j*(Omega_m[m]-Omega_m[n-1])*GFniwn[m])
        ContinuedFractionCoeff[n] = GFniwn[n]
#
#    for n in range(NumFreq):
#        print("%d      %12.10f   %12.10f")%(n, ContinuedFractionCoeff[n].real, ContinuedFractionCoeff[n].imag )     
#   
#   Rewrite C_M(z) = P_L(z)/Q_M(z) where 
#   P_L(z) = a_0 + a_1 z+ ...+ a_L z
#   Q_M(z) = b_0 + b_1 z+ ...+ b_M z  
#   with a_0=0, a_1=1, b_0=b_1=1, a_2=a_1, b_2=b_1+(z-z_1)*c_2
#   a_{i+1}(z) = a_i(z) + (z-z_i)*c_{i+1}*a_{i-1}(z)
#   b_{i+1}(z) = b_i(z) + (z-z_i)*c_{i+1}*b_{i-1}(z)  
#   Keep in mind that the bigger you chose N, i.e, the more coefficients c_i you calculate,
#   the bigger the numbers a_N and b_N are going to be.
#
    eta = 1.0e-8
    GFnw = numpy.zeros((len(EnergyMesh)), dtype=numpy.complex)
#   first algorithem to deal with large numbers
#    for i in range(len(EnergyMesh)):
#        pn1 = ContinuedFractionCoeff[0]
#        pn2 = ContinuedFractionCoeff[0]
#        pd1 = 1.0                                                                     
#        pd2 = 1.0+ContinuedFractionCoeff[1]*(EnergyMesh[i]+1j*eta-1j*Omega_m[0])     
#        for j in range(2,NumFreq):
#            jm = j-1
#            apj = ContinuedFractionCoeff[j]
#            cw = EnergyMesh[i]+1j*eta-1j*Omega_m[jm]
#            pn3 = pn2+apj*cw*pn1                 # P(energymesh,n)
#            pd3 = pd2+apj*cw*pd1                 # Q(energymesh,n)
#            confac = 1.0
#            if abs(pn3) > 1.0e+5: confac = 1.0*1.0e-5
#            pn1 = confac*pn2
#            pn2 = confac*pn3
#            pd1 = confac*pd2
#            pd2 = confac*pd3
#        GFnw[i] = pn3/pd3 
#                   
#    second algorithem to deal with large numbers 
    for i in range(len(EnergyMesh)):
        C1 = ContinuedFractionCoeff[0]          
        pn2 = 1.0
        pd2 = 1.0+ContinuedFractionCoeff[1]*(EnergyMesh[i]+1j*eta-1j*Omega_m[0])
        C2 = C1 * pn2 / pd2      
        for j in range(2,NumFreq):
            jm = j-1
            apj = ContinuedFractionCoeff[j]
            cw = EnergyMesh[i]+1j*eta-1j*Omega_m[jm]
            pn3 = 1.0+apj*cw/pn2                 # P(energymesh,n)
            pd3 = 1.0+apj*cw/pd2                 # Q(energymesh,n)
            C3 = C2 * pn3 / pd3            
#
            pn2 = pn3
            pd2 = pd3
            C2 = C3
        GFnw[i] = C3 
#
    return GFnw

if __name__ == "__main__":
    """ Analytic contniuation using Pade  """

    import argparse
#   setting up a parser
    parser = argparse.ArgumentParser()
#   When parse_args() is called, optional arguments will be identified by the - prefix, 
#   and the remaining arguments will be assumed to be positional
#   positional arguments
    parser.add_argument('filename', help=' filename ')
#   optional arguments
    parser.add_argument('-n', type=int, default=2, help=' column number ')
    parser.add_argument('-nfreq', type=int, help=' number of Freq. ')    
    parser.add_argument('-nene', type=int, default=200, help=' number of energy mesh ')
    parser.add_argument('-enecutoff', type=float, default = 2.0, help=' energy Cutoff ')
    parser.add_argument('-o', default='PadeOut.dat', help=' output filename ')   
#
    args = parser.parse_args()
    InFileName = args.filename
    ColumnNumber = args.n
    MaximumNumFreq = args.nfreq
    NumEnergyMesh = args.nene
    EnergyCutoff = args.enecutoff
    OutFileName = args.o    
    print "\n Input File Name: %s "%InFileName
    print " Column Number: %d "%ColumnNumber
    print " Number of energy mesh: %d "%NumEnergyMesh
    print " Energy cutoff: %12.6f "%EnergyCutoff
    print "\n Output File Name: %s "%OutFileName
#
    DataIn = ReadDataFile(InFileName, ColumnNumber)
    if not args.nfreq:
        MaximumNumFreq = len(DataIn[:,0])-1
    print " Num Freq.: %d "%MaximumNumFreq
    Omega_m = numpy.zeros((MaximumNumFreq), dtype=numpy.float)
    GFniwn = numpy.zeros((MaximumNumFreq), dtype=numpy.complex)
    Omega_m[0:MaximumNumFreq] = DataIn[0:MaximumNumFreq,0]
    GFniwn[0:MaximumNumFreq] = DataIn[0:MaximumNumFreq,1] + DataIn[0:MaximumNumFreq,2]*1j         
#
#   
    print "\n Input data:"
    for n in range(MaximumNumFreq):
        print("%12.10f      %12.10f   %12.10f")%(Omega_m[n], GFniwn[n].real, GFniwn[n].imag )
    if abs(GFniwn[0]) < 1.0e-10:
        print "\n Function is odd F(i\omega_0) = 0, the current implementation of continued fraction expansion is ill-defined"
        print " Comment out the first frequency in the input file and re-run the program!"
        quit()   
#
#   create energy mesh
    EnergyMesh = numpy.zeros((NumEnergyMesh+1), dtype=numpy.float)
    for n in range(NumEnergyMesh+1):
        EnergyMesh[n] = -EnergyCutoff + 2.0*float(n)*EnergyCutoff/float(NumEnergyMesh)    
#
    GFnw =  Pade(EnergyMesh, Omega_m, GFniwn)
#
    FileHandle = open(OutFileName, "w")    
    for n in range(NumEnergyMesh+1):
        FileHandle.write("%14.10f    %14.10f  %14.10f\n"%(EnergyMesh[n], GFnw[n].real, GFnw[n].imag ))
    FileHandle.close()
#
