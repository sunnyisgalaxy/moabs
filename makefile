##The static binary files bin/mcall, bin/mcomp, bin/numCI, bin/bsmap, should be able to run on ANY x86_64 Linux systems.
##"make install" command will compile the source codes on your system and will cp the dynamically linked binaries into bin/ to overwrite the static binary files. 
##You are suggested to compile the source code to use the version of BOOST, gcc, and other libs on your system.

MOABS_ROOT = $(shell pwd)

##ifndef BOOST_ROOT
##	BOOST_ROOT=$(MOABS_ROOT)/src/bsmap/boost
##endif

##ifndef SAMTOOLS
##	SAMTOOLS=$(MOABS_ROOT)/src/bsmap/samtools
##endif

CPAN=./.cpanm

all:
	@make -C src MOABS_ROOT=$(MOABS_ROOT) OPT=1
	$(CPAN) -l PERLLIBS Config::Simple threads threads::shared

install:
	@make -C src MOABS_ROOT=$(MOABS_ROOT) OPT=1 install

clean:
	@make -C src MOABS_ROOT=$(MOABS_ROOT) clean
.PHONY: clean

distclean: clean
	@rm -rf $(MOABS_ROOT)/bin
.PHONY: distclean
