FC=ifort -stand f08 -llapack -lblas -p -g -traceback -u\
	 -gen-interfaces -check all -warn declarations\
	 -warn interfaces -warn usage -fp-model strict
FC=ifort -stand f08 -llapack -lblas -p -g -traceback -u\
	 -gen-interfaces -check all -warn all -fp-model strict
FC77=ifort -llapack -lblas -O

ex1: movcol
	$(FC) -warn nounused -o $@ ex1.f90 libmovcol.a linpack/liblinpack.a

movcol: ddassl.o daux.o dlinpk.o movcol.o
	ar rs libmovcol.a movcol.o ddassl.o daux.o dlinpk.o linpack/liblinpack.a

movcol.o: movcol.f90
# at this point I don't want to remove the implicit declarations, so I
# turn the warnings off
	$(FC) -c $<

ddassl.o: ddassl.f90
	$(FC) -c $<

%.o: %.f
	$(FC77) -c $<

clean:
	rm -f *.o *.a
