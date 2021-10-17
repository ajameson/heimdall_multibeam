include $(top_srcdir)/cuda.mk

INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 

bin_PROGRAMS = \
	candidate_profiler \
	coincidencer \
	coincidencer_client \
	fil2pgm \
	generate_dmlist \
  heimdall

# To find the dedisp lib from bbarsdell
AM_CPPFLAGS = \
  -I$(top_srcdir)/Formats \
  -I$(top_srcdir)/Network \
  -I$(top_srcdir)/Pipeline \
  @PSRDADA_CFLAGS@ @DEDISP_CFLAGS@

heimdall_SOURCES = heimdall.C
coincidencer_SOURCES = coincidencer.C Candidates.C
coincidencer_client_SOURCES = coincidencer_client.C

LDADD = \
  $(top_builddir)/Formats/libhdformats.la \
  $(top_builddir)/Pipeline/libhdpipeline.la \
  $(top_builddir)/Network/libhdnetwork.la \
  @CUDA_LIBS@ @DEDISP_LIBS@ @PSRDADA_LIBS@

heimdall_CFLAGS   = $(CUDA_CFLAGS)
generate_dmlist_CXXFLAGS = @DEDISP_CFLAGS@

include $(top_srcdir)/config/Makefile.targets
 
