include $(top_srcdir)/cuda.mk

INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 

bin_PROGRAMS = heimdall coincidencer coincidencer_client candidate_profiler fil2pgm

# To find the dedisp lib from bbarsdell
AM_CPPFLAGS = -I../ -I../Formats -I../Network -I../Pipeline -I$(HOME) @PSRDADA_CFLAGS@ @DEDISP_CFLAGS@

heimdall_SOURCES = heimdall.C
coincidencer_SOURCES = coincidencer.C Candidates.C
coincidencer_client_SOURCES = coincidencer_client.C

LDADD = $(top_builddir)/Formats/libhdformats.la \
  $(top_builddir)/Pipeline/libhdpipeline.la \
  $(top_builddir)/Network/libhdnetwork.la \
  @CUDA_LIBS@ @DEDISP_LIBS@ @PSRDADA_LIBS@

heimdall_CFLAGS   = $(CUDA_CFLAGS)
#heimdall_LDADD    = $(CUDA_LIBS) @DEDISP_LIBS@ @PSRDADA_LIBS@
#heimdall_LDFLAGS  = $(CUDA_LDFLAGS)

include $(top_srcdir)/config/Makefile.targets

