bin_PROGRAMS = bseqc2

samtools_INCLUDE = $(top_srcdir)/lib/samtools
samtools_LIB = $(top_srcdir)/lib/samtools

CXXFLAGS = -g -O3 -std=c++11
bseqc2_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
bseqc2_LDFLAGS = -L$(samtools_LIB)
bseqc2_LDADD = -lbam -lz -lpthread -lboost_program_options -lboost_filesystem -lboost_system
bseqc2_SOURCES = bseqc2.cpp

bseqc2methplotdir = $(bindir)
bseqc2methplot_DATA = bseqc2mbiasplot.R
EXTRA_DIST = bseqc2mbiasplot.R
