# debug version
FCFLAGS=-r8 -i4 -g -p -fpe0 -traceback -llapack -lblas\
             -fp-model strict -ftrapuv -double-size 64
FC=ifort $(FCFLAGS) -u -gen-interfaces\
           -check all -warn all -warn nounused
FC77=ifort $(FCFLAGS)

# optimal version
# FCFLAGS=-r8 -i4 -O3 -traceback -llapack -lblas\
# #             -fp-model strict -ftrapuv
# FC=ifort $(FCFLAGS) -u -gen-interfaces\
# #           -warn all -warn nounused
# FC77=ifort $(FCFLAGS)

# gfortran
# FC=gfortran -llapack -lblas
# FC77=gfortran -llapack -lblas
