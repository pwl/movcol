libmovcol.so: movcol.o ddassl.o liblinpack.so
	@echo "Creating $@"
	@$(FC) -shared -o libmovcol.so movcol.o ddassl.o $(movcoldir)/linpack/liblinpack.so

movcol.o: movcol.f90 ddassl.o
	@echo "Compiling $@"
	@$(FC) -c $<

ddassl.o: ddassl.f90
	@echo "Compiling $@"
	@$(FC) -c $<

liblinpack.so:
	@$(MAKE) -C linpack

clean:
	@$(MAKE) clean -C linpack
	@rm -f *.o *.so *.mod
