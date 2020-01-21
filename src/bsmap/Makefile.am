bin_PROGRAMS = bsmap

gzstream_INCLUDE = $(top_srcdir)/lib/gzstream  
samtools_INCLUDE = $(top_srcdir)/lib/samtools

gzstream_LIB = $(top_srcdir)/lib/gzstream
samtools_LIB = $(top_srcdir)/lib/samtools

bsmap_CPPFLAGS = -DMAXHITS=1000 -DTHREAD -funroll-loops -I$(gzstream_INCLUDE) -I$(samtools_INCLUDE)
CXXFLAGS = -O3 -m64
bsmap_LDFLAGS = -L$(gzstream_LIB) -L$(samtools_LIB)
bsmap_LDADD = -lgzstream -lbam -lz -lpthread

bsmap_SOURCES =  \
	 param.h param.cpp \
	 utilities.h utilities.cpp \
	 dbseq.h dbseq.cpp \
	 reads.h reads.cpp \
	 align.h align.cpp \
	 pairs.h pairs.cpp \
	 main.cpp
