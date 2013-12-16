include makefile.inc

ex1: movcol
	$(FC) -o $@ ex1.f90 libmovcol.a linpack/liblinpack.a

harmonic: movcol
	$(FC) -o $@ harmonic.f90 libmovcol.a linpack/liblinpack.a

movcol: ddassl.o daux.o dlinpk.o movcol.o
	ar rs libmovcol.a movcol.o ddassl.o daux.o dlinpk.o linpack/liblinpack.a

movcol.o: movcol.f90
	$(FC) -c $<

ddassl.o: ddassl.f90
	$(FC) -c $<

%.o: %.f
	$(FC77) -c $<

clean:
	rm -f *.o *.a
