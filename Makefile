bindir = ../bin/
srcdir = ./src
# last edited 06/17/2018  M.S.Church
FC = mpif90 -O2
LIB = -llapack -lblas
OPTION = -mcmodel=large
.SUFFIXES: .f90 .o

#--------------------------------------------------
MOD = params.f90 \
      MonteCarlo.f90 \
      traj.f90 \
      supply.f90 
       

SRC = timecorr.f90

OBJ = $(MOD) $(SRC)

.f90.o:
	$(FC) -c $<
.f.o:
	$(FC) -c $<

dyn.x  :  $(OBJ)
	$(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

clear   :
	rm *.o
	rm *.x
	rm *.out 
	rm *.mod
	rm fort.*
