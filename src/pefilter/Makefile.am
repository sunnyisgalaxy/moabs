bin_PROGRAMS = pefilter

samtools_INCLUDE = $(top_srcdir)/lib/samtools
samtools_LIB = $(top_srcdir)/lib/samtools

CXXFLAGS = -g -O3 -std=c++11
pefilter_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
pefilter_LDFLAGS = -L$(samtools_LIB)
pefilter_LDADD = -lbam -lz -lpthread -lboost_program_options
pefilter_SOURCES = pefilter.cpp
