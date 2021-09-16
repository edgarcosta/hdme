# Move things to ./configure later on.

# Directories

# MODIFY THIS:
# CURDIR is where the source files lie
# MYDIR is where the other libraries are installed. See also ARB_DIR, etc.

CURDIR = /home/jean/code/modeq2
MYDIR = /home/jean/install

PARITWINE_DIR = $(MYDIR)
PARITWINE_INCDIR = $(PARITWINE_DIR)/include
PARITWINE_LIBDIR = $(PARITWINE_DIR)/lib

ARB_DIR = $(MYDIR)
ARB_INCDIR = $(ARB_DIR)/include
ARB_LIBDIR = $(ARB_DIR)/lib

FLINT_DIR = $(MYDIR)
FLINT_INCDIR = $(FLINT_DIR)/include/flint
FLINT_LIBDIR = $(FLINT_DIR)/lib

PARI_DIR = $(MYDIR)
PARI_INCDIR = $(PARI_DIR)/include/pari
PARI_LIBDIR = $(PARI_DIR)/lib

MPFR_DIR = $(MYDIR)
MPFR_INCDIR = $(MPFR_DIR)/include
MPFR_LIBDIR = $(MPFR_DIR)/lib

GMP_DIR = $(MYDIR)
GMP_INCDIR = $(GMP_DIR)/include
GMP_LIBDIR = $(GMP_DIR)/lib

# DO NOT MODIFY
# Compiling and linking

SHELL = /bin/sh

CC = gcc
CFLAGS = -ansi -Wall -pedantic -Werror -g

AR = ar
ARFLAGS = -rcs

VALGRIND = valgrind

INCS = -I$(CURDIR) -I$(PARITWINE_INCDIR) -I$(ARB_INCDIR) -I$(FLINT_INCDIR) \
-I$(PARI_INCDIR) -I$(MPFR_INCDIR) -I$(GMP_INCDIR)

LLIBS = -L$(PARITWINE_LIBDIR) -L$(ARB_LIBDIR) -L$(FLINT_LIBDIR) -L$(PARI_LIBDIR) \
-L$(MPFR_INCDIR) -L$(GMP_LIBDIR)

lDEPS = arb flint pari mpfr gmp pthread m # paritwine (at the beginning)
lLIBS = $(patsubst %, -l%, $(lDEPS))


# Modules to be compiled

MODULES = siegel hilbert pari isogeny acb_mat_extras theta igusa modular
SOURCES = $(wildcard $(patsubst %, %/*.c, $(MODULES)))
HEADERS = $(patsubst %, %.h, $(MODULES))
OBJS = $(patsubst %.c, build/%.o, $(SOURCES))

EXMP_SOURCES = $(wildcard examples/*.c)
EXMP_NAMES = $(patsubst %.c, %, $(EXMP_SOURCES))

TEST_SOURCES = $(wildcard $(patsubst %, %/test/t-*.c, $(MODULES)))
TESTS = $(patsubst %.c, build/test/%, $(TEST_SOURCES))
RUNTESTS = $(patsubst build/test/%, run/%, $(TESTS))
RUNVALGRIND = $(patsubst build/test/%, valgrind/%, $(TESTS))

TIMING_SOURCES = $(wildcard $(patsubst %, %/time/time-*.c, $(MODULES)))
TIMINGS = $(patsubst %.c, build/time/%, $(TIMING_SOURCES))

# Timing data is written in TIMEDIR

TIMEDIR = $(CURDIR)/time

# Targets

all: library

check: library tests $(RUNTESTS)

valgrind: library tests $(RUNVALGRIND)

tests: library $(TESTS)

timings: library $(TIMINGS)

clean:
	rm -rf build/
	rm libsea2.a

library: libsea2.a

libsea2.a: $(OBJS) $(SOURCES) $(HEADERS) | build
	$(AR) $(ARFLAGS) libsea2.a $(OBJS)

examples: library

build:
	@mkdir -p build build/test build/time \
	$(foreach mod, $(MODULES), build/$(mod)) \
	$(foreach mod, $(MODULES), build/test/$(mod)/test) \
	$(foreach mod, $(MODULES), build/time/$(mod)/time) 

build/%.o: %.c $(HEADERS) | build
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

run/%: build/test/%
	@$<

valgrind/%: build/test/%
	$(VALGRIND) --leak-check=full --track-origins=yes $<

build/test/%: %.c $(SOURCES) $(HEADERS) | build
	$(CC) $(CFLAGS) $(INCS) $(LLIBS) $< $(OBJS) -o $@ $(lLIBS)

build/time/%: %.c $(SOURCES) $(HEADERS) | build
	$(CC) -DTIMEDIR=\"$(TIMEDIR)\" $(CFLAGS) $(INCS) $(LLIBS) $< $(OBJS) -o $@ $(lLIBS) 

print-%:
	@echo '$*=$($*)'

plots:
	python $(TIMEDIR)/make_plots.py

gp2c:
	sage $(CURDIR)/gp2c/gp_poly_to_acb.sage

.PHONY: all check clean tests library examples build valgrind gp2c

.PRECIOUS: Makefile
