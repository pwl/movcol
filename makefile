MAKE  := make --no-print-directory
FCBASE:= ifort @$(PWD)/config/ifort.conf
FC    := $(FCBASE) @$(PWD)/config/ifort08.conf
FC77  := $(FCBASE) @$(PWD)/config/ifort77.conf

movcoldir := $(PWD)/movcol

export

all: movcol

examples/%: movcol
	@$(MAKE) $* -B -C examples

movcol:
	@$(MAKE) -C movcol

clean:
	@echo "Cleaning up"
	@$(MAKE) clean -C examples
	@$(MAKE) clean -C movcol

.PHONY: clean movcol
