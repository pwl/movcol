include makefile.inc

# ex1: movcol
# 	$(FC) -o $@ ex1.f90 libmovcol.a linpack/liblinpack.a

harmonic: movcol
	@echo "Compiling and linking harmonic"
	@$(FC) -o $@ harmonic.f90 libmovcol.a linpack/liblinpack.a

movcol: ddassl.o movcol.o
	@echo "Linking movcol"
	@ar rs libmovcol.a movcol.o ddassl.o linpack/liblinpack.a

movcol.o: movcol.f90
	@echo "Compiling movcol"
	@$(FC) -c $<

ddassl.o: ddassl.f90 linpack.a
	@echo "Compiling ddassl"
	@$(FC) -c $<

linpack.a:
	@echo "Compiling linpack"
	@$(MAKE) -C linpack

%.o: %.f
	$(FC77) -c $<

clean:
	rm -f *.o *.a
	$(MAKE) clean -C linpack
