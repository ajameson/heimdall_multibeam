
include $(top_srcdir)/cuda.mk

lib_LTLIBRARIES = libhdpipeline.la

libhdpipeline_la_SOURCES = \
  clean_filterbank_rfi.cu \
  default_params.C \
  error.C \
  find_giants.cu \
  get_rms.cu \
  label_candidate_clusters.cu \
  matched_filter.cu \
  measure_bandpass.cu \
  median_filter.cu \
  merge_candidates.cu \
  parse_command_line.C \
  pipeline.cu \
  remove_baseline.cu

nobase_include_HEADERS = \
  hd/error.h \
  hd/median_filter.h \
  hd/types.h

# To find the dedisp lib from bbarsdell
AM_CPPFLAGS = \
  -I$(top_srcdir)/Network \
  -I$(top_srcdir)/Formats \
  @DEDISP_CFLAGS@ @THRUST_CFLAGS@ @CUDA_CFLAGS@ $(BOOST_CPPFLAGS)

if HAVE_PSRDADA
AM_CPPFLAGS += @PSRDADA_CFLAGS@
endif

include $(top_srcdir)/config/Makefile.targets

