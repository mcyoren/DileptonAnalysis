
# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

# Fix missing runtime paths from RIVET.
ifeq ($(RIVET_USE),true)
  COMMA=,
  RIVET_LPATH=$(filter -L%,$(shell $(RIVET_BIN)$(RIVET_CONFIG) --ldflags))
  RIVET_RPATH=$(subst -L,-Wl$(COMMA)-rpath$(COMMA),$(RIVET_LPATH))
endif

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	$(info Usage: make mainXX)

# PYTHIA library.
$(PYTHIA):
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

WriteTreeOutput.so: WriteTreeOutputDct.cc
	$(CXX) $< -o $@ -w $(CXX_SHARED) $(CXX_COMMON)\
	 $(ROOT_LIB) `$(ROOT_CONFIG) --cflags`
WriteTreeOutputDct.cc: MyEvent.h WriteTreeOutputLinkDef.h
ifeq ($(ROOT_USE),true)
	$(ROOT_BIN)rootcint -f $@ -I$(PREFIX_INCLUDE) $^
else
	$(error Error: $@ requires ROOT)
endif

WriteTreeOutput: $(PYTHIA) $$@.cc WriteTreeOutput.so
	$(CXX) $@.cc WriteTreeOutput.so -o $@ -w $(CXX_COMMON) -Wl,-rpath,./\
	 `$(ROOT_CONFIG) --cflags --glibs`

WriteTreeoutputccbar: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
         `$(ROOT_CONFIG) --cflags --glibs`
else
        $(error Error: $@ requires ROOT)
endif

WriteTreeoutputbbbar WriteTreeOutputJetPairs WriteTreeoutputDrellYan: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
	`$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires ROOT)
endif

physics_processes_ccbar physics_processes_bbbar: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
	`$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires ROOT)
endif

Drell_Yan: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
         `$(ROOT_CONFIG) --cflags --glibs`
else
        $(error Error: $@ requires ROOT)
endif

main91 main91_pp main113_CuAu main113_pAu main113_pA main113_pA_nPDF main113_pp_nPDF: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
	 `$(ROOT_CONFIG) --cflags --glibs`\
      $(LHAPDF6_INCLUDE) $(LHAPDF6_LIB)
else
	$(error Error: $@ requires ROOT)
endif

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;\
	rm -f *.log; rm -f *plot.py; rm -f *.pcm; rm -f *.mpi;\
	rm -f *.hepmc; rm -f *.yoda; rm -f *.root;
