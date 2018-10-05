###############################################################################
#                                                                             #
#  make              ... generate executable                                  #
#  make clean        ... delete unnecessary files                             #
#                                                                             #
###############################################################################
#
#HOME    = ./
#SRCdir  = $(HOME)/SRC
#------------------------------------------------------------------------------
#
#
# FC ......... compiler name
# FOPT ....... optimization flags
# LDFLAGS .... linker flags
# LIBS ....... libraries needed to build the executable
# DESTDIR .... directory where the executable should go (incl. trailing '/'!)
# EXECNAME ... name of the resulting executable (without prefixed path!)
#
#------------------------------------------------------------------------------
# 
FC       = mpif90
FOPT     = -O3 -fPIC -xHOST -no-prec-div -check bounds -traceback # do not use -iop flag
CFLAGS   = -O3 -fPIC -xHOST -no-prec-div -check bounds -traceback # do not use -iop flag
CPP      = cpp -P -C -traditional
# For a serial run remove -DMPI switch in following line
CPPFLAGS = -DMPI
LDFLAGS  =
# 
LIBS	 =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
DESTDIR  = ../
EXECNAME = Tran
#
###############################################################################
###############################################################################
#
FFLAGS = $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)
#
#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS =  PrecisionMod.o ConstantsMod.o ParallelMod.o LibMod.o  \
	InputMod.o KPointsMod.o HamiltonianMod.o CurrentVertexMod.o SelfEnMod.o \
	OpticMod.o Main.o 
#
#..............................................................................
#
#  Build executable
#
$(EXEC):	$(OBJS) 
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS)
#..............................................................................
#
#  remove object files
#
clean:	
	rm  -f $(OBJS) 
#..............................................................................
#
#  define rules for generating object files from source files
#
.SUFFIXES : .o .f .for .f90
.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

.for.o:
	$(FC) $(FFLAGS) -c $<


