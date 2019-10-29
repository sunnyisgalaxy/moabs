bin_PROGRAMS = mcall

samtools_INCLUDE = $(top_srcdir)/lib/samtools
samtools_LIB = $(top_srcdir)/lib/samtools

mcall_CPPFLAGS = -Wall -w -I$(samtools_INCLUDE)
CXXFLAGS = -g -O3
mcall_LDFLAGS = -L$(samtools_LIB)


if OS_LINUX
mcall_LDADD = -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -lboost_date_time -lbam -lz -lrt -lpthread
endif

if OS_DARWIN
mcall_LDADD = -lboost_program_options -lboost_system -lboost_filesystem -lboost_thread -lboost_date_time -lbam -lz -lpthread
endif

mcall_SOURCES =  mcall.cpp
