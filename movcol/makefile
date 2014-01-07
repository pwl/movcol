libmovcol.so: movcol.o ddassl.o liblinpack.so libblas.so
	@echo "Creating $@"
	@$(FC) -shared -o libmovcol.so movcol.o ddassl.o $(movcoldir)/linpack/liblinpack.so $(movcoldir)/blas/libblas.so

movcol.o: movcol.f90 ddassl.o
	@echo "Compiling $@"
	@$(FC) -c $<

ddassl.o: ddassl.f90
	@echo "Compiling $@"
	@$(FC) -c $<

liblinpack.so:
	@echo "Compiling $@"
	@$(MAKE) -C linpack

libblas.so:
	@echo "Compiling $@"
	@$(MAKE) -C blas

clean:
	@$(MAKE) clean -C linpack
	@$(MAKE) clean -C blas
	@rm -f *.o *.mod
