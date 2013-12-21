# Linpack makefile
# MM
include ../makefile.inc

SRC     = $(wildcard *.f)
OBJ     = $(patsubst %.f,%.o,$(SRC))
LIB	= liblinpack.a

%.o: %.f
	@$(FC77) -c $<

linpack: $(OBJ)
	@echo "Creating liblinpack.a"
	@ar rs $(LIB) $(OBJ)

all: linpack

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *__genmod.f90

veryclean: clean
	rm -f $(LIB)
