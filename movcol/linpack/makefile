# Linpack makefile
# MM

SRC     = $(wildcard *.f)
OBJ     = $(patsubst %.f,%.o,$(SRC))
LIB	= liblinpack.so

$(LIB): $(OBJ)
	@$(FC77) -shared -o $(LIB) $(OBJ)

%.o: %.f
	@$(FC77) -c $<

.PHONY: clean

clean:
	@rm -f *.o *.mod *__genmod.f90 $(LIB)
